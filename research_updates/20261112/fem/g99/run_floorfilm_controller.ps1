param(
  [ValidateSet('MeshOnly','SolveApproved')][string]$Mode='MeshOnly',
  [string]$CaseId='G99_A3_P80_t3_l10p5_th0_phi0_s',
  [Parameter(Mandatory=$true)][string]$ExpectedManifestSha,
  [string]$ExpectedWindowSha='',
  [string]$StageRoot='__STAGE_ROOT__',
  [string]$ComsolBin='__COMSOL_BIN__',
  [string]$AttemptId='initial',
  [switch]$ResourceSelfTest,
  [switch]$UnitSelfTest,
  [switch]$PlanSelfTest,
  [switch]$EvaluationSelfTest
)
$ErrorActionPreference='Stop'
$OutputEncoding=[Console]::OutputEncoding=[Text.UTF8Encoding]::new($false)
$Cores=32;$SoftGiB=112;$HardGiB=128;$FreeGiB=100;$RuntimeFreeGiB=64;$GuardSamplesRequired=2
$ownedIds=@{};$ownedRootPid=$null;$ownedRootCreated=$null;$peakRss=0L;$active=$null
$manifestHash=$null;$admissionHash=$null;$windowHash=$null;$admissionDeadlineUtc=$null;$plan=$null;$case=$null;$lock=$null;$statusPath=$null
$started=[DateTime]::UtcNow.ToString('o');$completed=$null;$comparison=$null;$meshPrecheckRejection=$null

function Has-Field($Object,[string]$Name){return $null -ne $Object -and $null -ne $Object.PSObject.Properties[$Name]}
function Need($Object,[string]$Name){if(!(Has-Field $Object $Name)){throw ('Missing field: '+$Name)};return $Object.$Name}
function True($Value){return $Value -is [bool] -and $Value}
function Number($Value){try{$n=[double]$Value}catch{return $false};return $null -ne $Value -and $Value -isnot [bool] -and ![double]::IsNaN($n) -and ![double]::IsInfinity($n)}
function Num($Object,[string]$Name){$v=Need $Object $Name;if(!(Number $v)){throw ('Missing/nonfinite numeric evidence: '+$Name)};return [double]$v}
function Hash([string]$Path){if(!(Test-Path -LiteralPath $Path -PathType Leaf)){throw ('Missing file: '+$Path)};return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()}
function Json([string]$Path){return [IO.File]::ReadAllText($Path,[Text.Encoding]::UTF8)|ConvertFrom-Json}
function Write-Json([string]$Path,$Value){$tmp=$Path+'.tmp';[IO.File]::WriteAllText($tmp,($Value|ConvertTo-Json -Depth 40),[Text.UTF8Encoding]::new($false));Move-Item -LiteralPath $tmp -Destination $Path -Force}
function New-Json([string]$Path,$Value){if(Test-Path -LiteralPath $Path){throw ('Refusing to overwrite artifact: '+$Path)};Write-Json $Path $Value}
function Assert-EqualNumber($Object,[string]$Name,[double]$Expected){if([Math]::Abs((Num $Object $Name)-$Expected) -gt 1e-12){throw ('Fixed setting mismatch: '+$Name)}}
function Save-Status([string]$State,$Resources=$null,[string]$Reason=$null){
  if($null -eq $statusPath){return}

  if($null -ne $Resources){
    try{
      $now=[DateTime]::UtcNow
      if($null -eq $script:resourceJournalLastUtc){$script:resourceJournalLastUtc=[DateTime]::MinValue}
      if(($now-$script:resourceJournalLastUtc).TotalSeconds -ge 60 -or $State -like 'RESOURCE_STOP_*'){
        $row=[ordered]@{utc=$now.ToString('o');state=$State;case_id=$CaseId;mode=$Mode;free_GiB=$Resources.free_gib;host_CPU_percent=$Resources.cpu_percent;owned_RSS_GiB=$active.owned_rss_bytes/1GB;owned_process_ids=$active.owned_process_ids;host_low_samples=$active.host_low_samples;hard_samples=$active.hard_samples;D_free_GiB=$Resources.D_free_gib;D_healthy=$Resources.D_healthy;reason=$Reason}
        [IO.File]::AppendAllText((Join-Path $stage ('resources_'+$CaseId+'_'+$Mode+'.jsonl')),($row|ConvertTo-Json -Depth 8 -Compress)+[Environment]::NewLine,[Text.UTF8Encoding]::new($false));$script:resourceJournalLastUtc=$now
      }
      if($State -like 'RESOURCE_STOP_*'){
        $snap=[ordered]@{utc=$now.ToString('o');state=$State;case_id=$CaseId;mode=$Mode;resources=$Resources;active=$active;peak_owned_rss_bytes=$peakRss;processes=@(Get-CimInstance Win32_Process|Select-Object ProcessId,ParentProcessId,Name,WorkingSetSize,CreationDate)}
        New-Json (Join-Path $stage ('guard_stop_'+$CaseId+'_'+$now.ToString('yyyyMMddTHHmmssfff')+'.json')) $snap
      }
    }catch{[Console]::Error.WriteLine('RESOURCE_JOURNAL_WARNING: '+$_.Exception.Message)}
  }

  Write-Json $statusPath ([ordered]@{schema=1;state=$State;case_id=$CaseId;mode=$Mode;attempt_id=$AttemptId;started_utc=$started;updated_utc=[DateTime]::UtcNow.ToString('o');active=$active;completed=$completed;comparison=$comparison;reason=$Reason;mesh_precheck_rejection=$meshPrecheckRejection;resource=$Resources;cores=$Cores;concurrency=1;peak_owned_rss_bytes=$peakRss;owned_RSS_soft_GiB=$SoftGiB;owned_RSS_hard_GiB=$HardGiB;startup_free_GiB=$FreeGiB;runtime_free_GiB=$RuntimeFreeGiB;guard_consecutive_samples=$GuardSamplesRequired;manifest_sha256=$manifestHash;admission_sha256=$admissionHash;queue_window_sha256=$windowHash;new_admission_deadline_utc=$admissionDeadlineUtc;automatic_retries=0;solver_wall_timeout_s=$null;cumulative_budget_s=$null;scope='G99 P80 10.5um normal-s 3um bottom fill. Raw net flux warnings remain separate from accepted points.'})
}
function Assert-Manifest {
  $path=Join-Path $stage 'launch_manifest.json';$h=Hash $path
  if($h -ne $ExpectedManifestSha.ToLowerInvariant() -or ($null -ne $script:manifestHash -and $h -ne $script:manifestHash)){throw 'Manifest differs from explicit expected/frozen SHA256'}
  $m=Json $path;$sources=Need $m 'source_hashes'
  foreach($n in @('CavityFloorFilm.java','run_floorfilm_controller.ps1','plan.json','material.json','warn_fixture_metrics.json')){[void](Need $sources $n)}
  foreach($p in $sources.PSObject.Properties){$n=[string]$p.Name;$v=[string]$p.Value;if([IO.Path]::GetFileName($n) -ne $n -or $n -match '[\\/]' -or $n -eq 'launch_manifest.json' -or $v -notmatch '^[a-f0-9]{64}$' -or (Hash (Join-Path $stage $n)) -ne $v){throw ('Source fingerprint mismatch/unsafe name: '+$n)}}
  if([IO.Path]::GetFullPath($PSCommandPath) -ne (Join-Path $stage 'run_floorfilm_controller.ps1')){throw 'Controller must execute its manifest-bound stage copy'}
  $script:manifestHash=$h
}
function Case-Contract([string]$Id){
    $index = [Array]::IndexOf(@('G99_A3_P80_t3_l10p5_th0_phi0_s','G99_A0facets_P80_t3_l10p5_th0_phi0_s'), $Id)
    if ($index -lt 0) { throw 'Only the two approved G99 cases are executable' }
    return [pscustomobject]@{
        height_um=98.09681325894577; opening_flat_to_flat_um=71.74019438588937; bottom_flat_to_flat_um=34.65468216965557; bottom_width_um=34.65468216965557
        film_um=3; lambda_um=10.5; theta_deg=0; phi_deg=0; polarization='s'; mesh='m4_air29_wft3'
        source_geometry_id=@('G99_A3','G99_A0_facets')[$index]
        wall_profile_id=@('A3_SINGLE_STRAIGHT_PANEL_PER_SIDE','A0_EIGHT_PLANAR_FACETS_PER_SIDE')[$index]
        geometry_marker=@('G99_A3_STRAIGHT_P80_SIO2_FLOOR_3UM','G99_A0_8FACET_P80_SIO2_FLOOR_3UM')[$index]
        profile_segments=@(1,8)[$index]
    }
}

function Source-Phase([string]$Id){[void](Case-Contract $Id);return 'Main'}
function Flat-Path([string]$Id){[void](Case-Contract $Id);return (Join-Path $stage 'flat_control_main_admission.json')}
function Flat-State([string]$Id){[void](Case-Contract $Id);return 'PASS_G99_MAIN_FLAT_SOURCE_ADMISSION'}
function Flat-Count([string]$Id){[void](Case-Contract $Id);return 2}
function Mesh-FamilyKey($Row){return ((Num $Row 'height_um').ToString('R',[Globalization.CultureInfo]::InvariantCulture)+'|'+(Need $Row 'mesh')+'|'+(Num $Row 'bottom_flat_to_flat_um').ToString('R',[Globalization.CultureInfo]::InvariantCulture)+'|'+(Num $Row 'film_um').ToString('R',[Globalization.CultureInfo]::InvariantCulture))}
function Case-Polarization([string]$Id){return (Case-Contract $Id).polarization}

function Tangential-ModeAxis($Row){
 $p=(Num $Row 'phi_deg')*[Math]::PI/180.;$pol=[string](Need $Row 'polarization')
 if($pol -eq 'p'){return @([Math]::Cos($p),[Math]::Sin($p),0.)}
 if($pol -eq 's'){return @((-[Math]::Sin($p)),[Math]::Cos($p),0.)}
 throw 'Only explicit p/s tangential source axes allowed'
}

function Case-FilmThickness([string]$Id){return (Case-Contract $Id).film_um}
function Incidence($Row){
 $theta=(Num $Row 'theta_deg')*[Math]::PI/180.;$phi=(Num $Row 'phi_deg')*[Math]::PI/180.;$k0=2*[Math]::PI/((Num $Row 'lambda_um')*1e-6)
 $st=[Math]::Sin($theta);$ct=[Math]::Cos($theta);$sp=[Math]::Sin($phi);$cp=[Math]::Cos($phi);$pol=[string](Need $Row 'polarization')
 if($pol -eq 'p'){$ein=@(($ct*$cp),($ct*$sp),$st);$eout=@(($ct*$cp),($ct*$sp),(-$st))}elseif($pol -eq 's'){$ein=@((-$sp),$cp,0.);$eout=@((-$sp),$cp,0.)}else{throw 'Only explicit p/s polarization allowed'}
 return [ordered]@{k0=$k0;sin_theta=$st;cos_theta=$ct;k_in=@(($k0*$st*$cp),($k0*$st*$sp),(-$k0*$ct));k_out=@(($k0*$st*$cp),($k0*$st*$sp),($k0*$ct));E_in=$ein;E_out=$eout}
}
function Assert-Vector($Actual,$Expected,[double]$Tolerance,[string]$Name){
  $a=@($Actual);$e=@($Expected);if($a.Count -ne 3 -or $e.Count -ne 3){throw ('Expected three-component vector: '+$Name)}
  for($j=0;$j -lt 3;$j++){if(!(Number $a[$j]) -or [Math]::Abs([double]$a[$j]-[double]$e[$j]) -gt $Tolerance){throw ('Vector evidence mismatch: '+$Name)}}
}
function Expected-Lattice($Row,[double]$FrequencyRatio=1.){
  $inc=Incidence $Row;$k0=$inc.k0*$FrequencyRatio;$bx=2*[Math]::PI/80e-6;$by=$bx/[Math]::Sqrt(3);$orders=@()
  for($a=-40;$a -le 40;$a++){for($b=-40;$b -le 40;$b++){$x=$inc.k_in[0]*$FrequencyRatio+$a*$bx;$y=$inc.k_in[1]*$FrequencyRatio+(-$a+2*$b)*$by;if(($x*$x+$y*$y)/($k0*$k0) -lt 1-1e-12){$orders+=,[pscustomobject]@{m=$a;n=$b;kx=$x;ky=$y}}}}
  return $orders
}
function Assert-LatticeIndices($Modes,[string]$Name,$Independent){
  $indices=@(Need $Modes $Name);$orders=@($Independent)
  if($indices.Count -ne $orders.Count){throw 'Frozen lattice order index count differs from independent incidence-shifted lattice'}
  for($j=0;$j -lt $indices.Count;$j++){$pair=@($indices[$j]);if($pair.Count -ne 2 -or !(Number $pair[0]) -or !(Number $pair[1]) -or $pair[0] -ne $orders[$j].m -or $pair[1] -ne $orders[$j].n){throw 'Frozen lattice order indices differ from independent incidence-shifted lattice'}}
}
function Assert-ExecutableCase([string]$Id){[void](Case-Contract $Id)}
function Assert-FitDirection($Metrics,$Fit,$Row){
  $inc=Incidence $Row;$k=@(Need $Metrics 'incident_k_rad_per_m');Assert-Vector $k $inc.k_in .1 'physical incoming k before normalized-direction power-fit consistency'
  $norm2=0.;foreach($v in $k){$norm2+=[double]$v*[double]$v};$norm=[Math]::Sqrt($norm2)
  if($norm -le 0 -or [Math]::Abs($norm-$inc.k0) -gt .1){throw 'Actual physical vacuum wavevector magnitude differs from requested frequency'}
  $ct=-[double]$k[2]/$norm;if($ct -le 0 -or $ct -gt 1){throw 'Actual incoming wavevector does not point downward'}
  if(!(True (Need $Fit 'direction_normalized_to_exact_vacuum_k0'))){throw 'Fit must normalize actual direction to the exact vacuum k0'}
  Assert-Vector (Need $Fit 'actual_k_in_rad_per_m') $k 1e-8 'fit actual raw incoming wavevector metadata'
  $projected=@();foreach($v in $k){$projected+=([double]$v*$inc.k0/$norm)};Assert-Vector (Need $Fit 'projection_k_in_rad_per_m') $projected 1e-6 'fit normalized exact-frequency wavevector metadata'
  if([Math]::Abs((Num $Fit 'actual_k_magnitude_rad_per_m')/$norm-1) -gt 1e-12){throw 'Fit wavevector-magnitude metadata differs from actual measured vector'}
  if([Math]::Abs((Num $Fit 'cos_theta')-$ct) -gt 1e-12 -or [Math]::Abs((Num $Fit 'k0_fixed_rad_per_m')-$inc.k0) -gt 1e-7 -or [Math]::Abs((Num $Fit 'kz_fixed_rad_per_m')-$inc.k0*$ct) -gt 1e-7 -or [Math]::Abs((Num $Fit 'power_admittance_p_over_vacuum')-1/$ct) -gt 1e-12 -or [Math]::Abs((Num $Fit 'power_admittance_s_over_vacuum')-$ct) -gt 1e-12){throw 'Fit direction/admittance metadata inconsistent with minus-kz-over-actual-vector-norm'}
  return $ct
}
function Assert-PresolvePortProof($Proof){
  if((Need $Proof 'case_id') -ne $CaseId -or (Need $Proof 'proof_schema') -ne 'ANGLE_BAND_SOURCE_POWER_V2' -or (Need $Proof 'state') -ne 'PASS_PRESOLVE_COMPLETE_PORT_BASIS' -or (Need $Proof 'solver_executed') -isnot [bool] -or $Proof.solver_executed){throw 'Actual case-bound Variables-stage complete-port proof required before physical solve'}
  foreach($n in @('physical_frequency_restored','incident_wavevector_matches_requested','actual_wavevectors_match_independent_lattice','physical_index_set_exact','padding_index_set_exact','port_names_unique')){if(!(True (Need $Proof $n))){throw ('Presolve port proof failed: '+$n)}}
  if((Need $Proof 'source_frequency_expression') -ne 'f0'){throw 'Presolve source frequency was not restored to the physical f0'}
  Assert-EqualNumber $Proof 'generator_superset_frequency_ratio' 1.2;Assert-EqualNumber $Proof 'retained_padding_frequency_ratio' 1.1
  $physical=@(Expected-Lattice $case);$padding=@(Expected-Lattice $case 1.1);$modes=Need $case 'expected_modes';$inc=Incidence $case;$k0=$inc.k0;$bx=2*[Math]::PI/80e-6;$by=$bx/[Math]::Sqrt(3)
  Assert-EqualNumber $Proof 'physical_spatial_orders' $physical.Count;Assert-EqualNumber $Proof 'physical_polarization_channels' (2*$physical.Count);Assert-EqualNumber $Proof 'total_port_channels' (2*$padding.Count);Assert-EqualNumber $Proof 'evanescent_channels' (2*($padding.Count-$physical.Count))
  $counts=@(Need $Proof 'actual_independent_spatial_channel_counts');if($counts.Count -ne $physical.Count -or @($counts|Where-Object {!(Number $_) -or $_ -ne 2}).Count){throw 'Presolve physical channels do not include exactly two polarizations per independent spatial order'}
  $expectedPadding=@{};$expectedPhysical=@{};foreach($p in $padding){$expectedPadding[([string]$p.m+','+[string]$p.n)]=0};foreach($p in $physical){$expectedPhysical[([string]$p.m+','+[string]$p.n)]=0}
  $rows=@(Need $Proof 'channels');if($rows.Count -ne 2*$padding.Count){throw 'Presolve total actual port inventory differs from frozen exact padding set'}
  $ports=@{};$keys=@{};$components=@{};$sourceCount=0
  foreach($r in $rows){
    $pn=[string](Need $r 'port_name');$key=[string](Need $r 'key');if([string]::IsNullOrWhiteSpace($pn) -or $ports.ContainsKey($pn) -or $keys.ContainsKey($key)){throw 'Presolve port names or channel keys duplicate'};$ports[$pn]=$true;$keys[$key]=$true
    $a=Num $r 'canonical_m';$b=Num $r 'canonical_n';if($a -ne [Math]::Floor($a) -or $b -ne [Math]::Floor($b)){throw 'Presolve canonical lattice indices must be integers'};$pair=[string]([int]$a)+','+[string]([int]$b)
    if(!$expectedPadding.ContainsKey($pair)){throw 'Presolve actual port outside the original exact 1.1 padding index set'};$expectedPadding[$pair]++
    if(!(True (Need $r 'wavevector_verified'))){throw 'Presolve actual global port wavevector not verified'}
    $kx=$inc.k_in[0]+$a*$bx;$ky=$inc.k_in[1]+(-$a+2*$b)*$by;$radial=$kx*$kx+$ky*$ky;$prop=$radial -lt $k0*$k0*(1-1e-12)
    foreach($field in @('propagating','actual_wavevector_propagating')){if((Need $r $field) -isnot [bool] -or $r.$field -ne $prop){throw 'Presolve physical/evanescent classification inconsistent with actual fixed-frequency lattice'}}
    if([Math]::Abs((Num $r 'kx_rad_per_m')-$kx) -gt .1 -or [Math]::Abs((Num $r 'ky_rad_per_m')-$ky) -gt .1){throw 'Presolve actual transverse wavevector differs from independent canonical reciprocal lattice'}
    $rz=Num $r 'real_kz_rad_per_m';$iz=Num $r 'imag_kz_rad_per_m'
    if($prop){if([Math]::Abs($iz) -gt .1 -or [Math]::Abs([Math]::Abs($rz)-[Math]::Sqrt($k0*$k0-$radial)) -gt .1 -or ($key -ne '0,0,co' -and $rz -le 0) -or !$expectedPhysical.ContainsKey($pair)){throw 'Presolve propagating kz or independent physical membership mismatch'};$expectedPhysical[$pair]++}
    else{if([Math]::Abs($rz) -gt .1 -or [Math]::Abs([Math]::Abs($iz)-[Math]::Sqrt($radial-$k0*$k0)) -gt .1){throw 'Presolve evanescent kz differs from fixed-frequency vacuum dispersion'}}
    if($key -notmatch '^(-?\d+),(-?\d+),(.+)$'){throw 'Presolve malformed native order/component key'};$component=$Matches[3];if(!$components.ContainsKey($pair)){$components[$pair]=@()};$components[$pair]+=$component
    if($key -eq '0,0,co'){if($pn -ne '1' -or $pair -ne '0,0' -or !$prop){throw 'Presolve source channel identity changed'};$sourceCount++}
  }
  if($sourceCount -ne 1){throw 'Exactly one actual source co channel required'}
  foreach($pair in $expectedPadding.Keys){if($expectedPadding[$pair] -ne 2){throw 'Presolve padding basis is incomplete or duplicates a spatial order'};$wanted=if($pair -eq '0,0'){'co|orthogonal'}else{'inplane|outofplane'};if((@($components[$pair]|Sort-Object) -join '|') -ne $wanted){throw 'Presolve exact spatial order lacks two distinct polarization channels'}}
  foreach($pair in $expectedPhysical.Keys){if($expectedPhysical[$pair] -ne 2){throw 'Presolve physical basis is incomplete or duplicates a spatial order'}}
}
function Case-GeometryMarker { return (Case-Contract $CaseId).geometry_marker }

function Analytic-Geometry($Row) {
    $z=@($Row.profile_knot_heights_um); $a=@($Row.profile_knot_apothems_um)
    $segments=[int](Case-Contract $Row.id).profile_segments
    if ($z.Count -ne $segments+1 -or $a.Count -ne $z.Count -or $z[0] -ne 0 -or
        $z[-1] -ne $Row.height_um -or [Math]::Abs(2*$a[0]-$Row.bottom_flat_to_flat_um) -gt 1e-10 -or
        [Math]::Abs(2*$a[-1]-$Row.opening_flat_to_flat_um) -gt 1e-10) {
        throw 'Frozen G99 profile does not bind the approved endpoints'
    }
    $t=[double]$Row.film_um
    if ($z[1] -le $t) { throw 'Film must cut the first linear panel' }
    $c=2*[Math]::Sqrt(3); $cell=$c*40*40; $floor=$c*$a[0]*$a[0]; $open=$c*$a[-1]*$a[-1]
    $fill=$a[0]+($a[1]-$a[0])*$t/$z[1]; $film=$c*$t*($a[0]*$a[0]+$a[0]*$fill+$fill*$fill)/3
    $total=$cell*10.5; $wall=0.
    for ($i=1; $i -lt $z.Count; $i++) {
        $dz=$z[$i]-$z[$i-1]; $da=$a[$i]-$a[$i-1]
        if ($dz -le 0 -or $da -le 0) { throw 'Non-monotone frozen G99 wall profile' }
        $total+=$c*$dz*($a[$i-1]*$a[$i-1]+$a[$i-1]*$a[$i]+$a[$i]*$a[$i])/3
        $wall+=$c*($a[$i-1]+$a[$i])*[Math]::Sqrt($dz*$dz+$da*$da)
    }
    return [ordered]@{
        floor_fill_thickness_um=$t; top_area_m2=$cell*1e-12; floor_area_m2=$floor*1e-12
        land_area_m2=($cell-$open)*1e-12; wall_area_m2=$wall*1e-12
        film_volume_m3=$film*1e-18; expected_film_volume_m3=$film*1e-18
        air_volume_m3=($total-$film)*1e-18; total_computational_volume_m3=$total*1e-18
        expected_total_volume_m3=$total*1e-18; film_top_apothem_um=$fill; film_air_interface_area_m2=$c*$fill*$fill*1e-12
    }
}

function Assert-ProfileBinding($Geometry,$Row) {
    $contract=Case-Contract $Row.id
    foreach ($name in @('source_geometry_id','wall_profile_id')) {
        if ((Need $Geometry $name) -cne (Need $contract $name) -or (Need $Row $name) -cne (Need $contract $name)) {
            throw ('Native G99 wall profile differs: '+$name)
        }
    }
    if ((Need $Geometry 'geometry_parameter_fingerprint_scheme') -cne 'SHA256_G99_HEX_PIECEWISE_FLOOR_FILL_V1_HEX_DOUBLE_SI' -or
        (Need $Geometry 'geometry_parameter_fingerprint_sha256') -cnotmatch '^[a-f0-9]{64}$') { throw 'G99 profile fingerprint missing' }
    foreach ($names in @(@('profile_knot_heights_m','profile_knot_heights_um'),@('profile_knot_apothems_m','profile_knot_apothems_um'))) {
        $actual=@(Need $Geometry $names[0]); $expected=@(Need $Row $names[1])
        if ($actual.Count -ne $expected.Count) { throw 'G99 profile knot count differs' }
        for ($j=0; $j -lt $actual.Count; $j++) {
            if (!(Number $actual[$j]) -or [Math]::Abs($actual[$j]-1e-6*$expected[$j]) -gt 1e-14) {
                throw ('Native G99 profile knot differs: '+$names[0])
            }
        }
    }
}

function Assert-PairAdmission {
    $path=Join-Path $stage 'pair_mesh_admission.json'
    $pair=Json $path
    if ($pair.state -ne 'PASS_G99_TWO_MESH_PREFLIGHTS' -or $pair.stage_root -ne $stage -or
        $pair.manifest_sha256 -ne $manifestHash -or $pair.queue_window_sha256 -ne $windowHash -or
        (@($pair.entries.case_id) -join '|') -ne ('G99_A3_P80_t3_l10p5_th0_phi0_s|G99_A0facets_P80_t3_l10p5_th0_phi0_s')) {
        throw 'Both G99 mesh preflights must pass before either physical solve'
    }
    foreach ($entry in $pair.entries) {
        $bound=Join-Path $stage ('admission_'+$entry.case_id+'.json')
        if ($entry.admission_sha256 -ne (Hash $bound)) { throw 'Pair mesh admission was altered' }
    }
}

function Admission-Path {return (Join-Path $stage ('admission_'+$CaseId+'.json'))}

function Initialize-ExecutionContext {
 Assert-ExecutableCase $CaseId
 Assert-Manifest
 Assert-Window
 $script:plan=Json (Join-Path $stage 'plan.json')
 Assert-Plan
 [void](Assert-ConditionalSourceAdmission)
}

function Assert-Window {
 if($ExpectedWindowSha -notmatch '^[a-fA-F0-9]{64}$'){throw 'Explicit execution-context SHA required'}
 $path=Join-Path $stage 'window.json';$h=Hash $path;if($h -ne $ExpectedWindowSha.ToLowerInvariant() -or ($null -ne $script:windowHash -and $h -ne $script:windowHash)){throw 'Immutable context changed'};$w=Json $path
 if($w.purpose -ne 'G99_PAIR_UNWINDOWED_CONTEXT' -or $w.stage_root -ne $stage -or $w.manifest_sha256 -ne $manifestHash -or $w.time_window_enforced -ne $false -or $null -ne $w.new_admission_deadline_utc -or $w.automatic_retries -ne 0){throw 'Closed unwindowed context binding failed'}
 $script:windowHash=$h;$script:admissionDeadlineUtc=$null
}
function Assert-OOC($Override){
  if((Need $Override 'case_id') -ne $CaseId -or (Need $Override 'state') -ne 'VERIFIED_OOC_CONFIGURATION_BEFORE_PHYSICS' -or (Need $Override 'method') -ne 'EXACT_DIRECT_OUT_OF_CORE_NO_BLR' -or (Need $Override 'solver_tag') -ne 'sol1' -or (Need $Override 'stationary_tag') -ne 's1'){throw 'Actual case-bound production OOC readback missing'}
  foreach($n in @('readback_verified','active_selection_preserved','oocmemory_is_not_process_RSS_cap')){if(!(True (Need $Override $n))){throw ('OOC invariant failed: '+$n)}}
  if((Need $Override 'physical_thresholds_changed') -isnot [bool] -or $Override.physical_thresholds_changed){throw 'OOC must not change physical thresholds'}
  foreach($n in @('active_selection_before','active_selection_after')){$v=Need $Override $n;if(@($v.PSObject.Properties).Count -ne 1 -or (Need $v 'fc1') -ne 'd1'){throw 'The actual active fc1->d1 direct solver selection must be preserved'}}
  $rows=@(Need $Override 'direct_changes');if($rows.Count -ne 2 -or (Num $Override 'direct_feature_count') -ne 2){throw 'Exactly both existing main Direct features require verified OOC'}
  $seen=@{}
  foreach($r in $rows){
    $tag=[string](Need $r 'tag');if($tag -notin @('dDef','d1') -or $seen.ContainsKey($tag) -or !(True (Need $r 'readback_verified'))){throw 'Duplicate/missing OOC Direct feature'};$seen[$tag]=$true
    $kind=if($tag -eq 'dDef'){'mumps'}else{'pardiso'};$before=Need $r 'before';$after=Need $r 'after'
    if((Need $before 'linsolver') -ne $kind -or (Need $after 'linsolver') -ne $kind -or (Need $after 'ooc') -ne 'on' -or (Need $after 'incore') -ne 'manual' -or (Need $after 'mumpsblr') -ne 'off'){throw 'Direct solver kind/OOC/manual/no-BLR contract changed'}
    Assert-EqualNumber $after 'oocmemory_MB' 24576
  }
}
function Assert-Plan {
  $g=Need $script:plan 'geometry'
  foreach($x in @(@('pitch_um',80),@('opening_flat_to_flat_um',71.74019438588937),@('air_cap_um',10.5))){Assert-EqualNumber $g $x[0] $x[1]}
  if((Need $g 'structure') -ne 'G99_straight_facets_matched_floorfill' -or (Need $g 'height_um') -ne 98.09681325894577 -or (Need $g 'film_um') -ne 3 -or (Need $g 'film_configuration') -ne 'floorfill_only_not_conformal' -or (Need $g 'polarization') -ne 'true_global_case_azimuth_p_or_s'){throw 'Only the matched G99 wall pair is authorized'}
  $tol=Need $script:plan 'acceptance'
  $fixed=@{absolute_power_mismatch_max=.005;relative_absorption_mismatch_max=.02;incident_power_relative_error_max=.005;incident_fit_relative_L2_residual_max=.005;incident_fit_leave_one_out_Pin_error_max=.005;side_pair_power_mismatch_over_Pin_max=.002;passivity_abs_tolerance=.002;strict_absolute_power_mismatch_max=.002;incident_fit_condition_number_max=10;actual_source_polarization_leakage_fraction_max=.000001}
  foreach($n in $fixed.Keys){Assert-EqualNumber $tol $n $fixed[$n]}
  $resources=Need $script:plan 'resource_plan'
  $limits=@{cores_per_case=32;concurrency=1;owned_RSS_soft_guard_GiB=112;owned_RSS_hard_guard_GiB=128;Windows_free_before_admission_GiB_min=100;host_free_memory_runtime_minimum_GiB=64;host_CPU_before_admission_percent_max=75;guard_consecutive_samples=2;actual_DOF_max=8000000;tetrahedra_max=2500000;min_tetrahedron_quality=.2;D_free_GiB_min=100}
  foreach($n in $limits.Keys){Assert-EqualNumber $resources $n $limits[$n]}
  $rows=@(Need $script:plan 'cases')
  $wanted=@('G99_A3_P80_t3_l10p5_th0_phi0_s','G99_A0facets_P80_t3_l10p5_th0_phi0_s')
  if($rows.Count -ne $wanted.Count -or (@($rows|ForEach-Object {[string](Need $_ 'id')}) -join '|') -ne ($wanted -join '|')){throw 'Exactly the two ordered G99 cases are authorized'}
  foreach($row in $rows){
    $contract=Case-Contract ([string](Need $row 'id'));if((Need $row 'polarization') -ne $contract.polarization -or (Need $row 'mesh') -ne $contract.mesh){throw 'Case id/polarization/mesh binding mismatch'}
    foreach($n in @('height_um','opening_flat_to_flat_um','bottom_flat_to_flat_um','film_um','lambda_um','theta_deg','phi_deg')){Assert-EqualNumber $row $n $contract.$n};$thickness=$contract.film_um
    foreach($n in @('n_film','k_film','n_ti','k_ti')){if((Num $row $n) -le 0){throw 'Passive source-proxy film/Ti n,k must be positive'}}
    $analytic=Analytic-Geometry $row;$expect=Need $row 'expected_geometry'
    $geometryMapping=@{film_top_apothem_um=@('film_top_apothem_um',1.);film_volume_um3=@('film_volume_m3',1e18);total_dielectric_volume_um3=@('total_computational_volume_m3',1e18);air_volume_um3=@('air_volume_m3',1e18);film_air_interface_area_um2=@('film_air_interface_area_m2',1e12);top_area_um2=@('top_area_m2',1e12)}
    foreach($n in $geometryMapping.Keys){$bound=$geometryMapping[$n];$want=[double]$analytic[$bound[0]]*[double]$bound[1];if([Math]::Abs((Num $expect $n)/$want-1) -gt 1e-10){throw ('Case expected geometry is not the independent analytic value: '+$n)}}
    foreach($x in @(@('material_domains',2),@('periodic_pairs',3),@('film_hmax_um',$(if($contract.lambda_um -eq 9){.2}else{.5})),@('air_hmax_um',2.9),@('Ti_hmax_um',$(if($contract.lambda_um -eq 9){.75}else{.875})),@('edge_hmax_um',$(if($contract.lambda_um -eq 9){.225}else{.2625})))){Assert-EqualNumber $expect $x[0] $x[1]}
    if($contract.mesh -eq 'm4_air29_wft3_localtop24'){Assert-EqualNumber $expect 'top_port_hmax_um' 2.4}elseif($null -ne (Need $expect 'top_port_hmax_um')){throw 'Unapproved local top mesh in coarse case'}
    $modes=Need $row 'expected_modes';Assert-EqualNumber $modes 'padding_generation_frequency_ratio' 1.10
    $physical=@(Expected-Lattice $row);$padding=@(Expected-Lattice $row 1.10)
    Assert-EqualNumber $modes 'spatial_propagating_orders' $physical.Count;Assert-EqualNumber $modes 'propagating_polarization_channels' (2*$physical.Count);Assert-EqualNumber $modes 'all_polarization_channels' (2*$padding.Count)
    Assert-LatticeIndices $modes 'physical_order_indices' $physical;Assert-LatticeIndices $modes 'padding_order_indices' $padding
    $planes=@(Need $expect 'fixed_uniform_air_source_plane_z_um');$height=Num $row 'height_um';$zWant=@(($height+1.05),($height+2.625),($height+5.25),($height+7.875),($height+9.45));if($planes.Count -ne 5){throw 'Exactly five fixed source planes required'};for($j=0;$j -lt 5;$j++){if(!(Number $planes[$j]) -or [Math]::Abs([double]$planes[$j]-$zWant[$j]) -gt 1e-10){throw 'Fixed source-plane inventory changed'}}
  }
  foreach($row in $rows){$same=@($rows|Where-Object {$_.lambda_um -eq $row.lambda_um});foreach($r in $same){foreach($n in @('n_film','k_film','n_ti','k_ti')){if([Math]::Abs((Num $r $n)-(Num $row $n)) -gt 1e-12){throw 'Same-wavelength thickness/angle controls must use exactly the same optical material'}}}}
  $selected=@($rows|Where-Object {$_.id -eq $CaseId});if($selected.Count -ne 1){throw 'Requested case outside the frozen case set'};$script:case=$selected[0]
  $mat=Json (Join-Path $stage 'material.json')
  if((Need $mat 'kind') -ne 'PUBLIC_OPTICAL_PROXY_CASE_SELECTED_FLOORFILL_SIO2_ON_TI' -or !(True (Need $mat 'case_selected'))){throw 'Case-selected traceable optical proxy material receipt required'}
  $optical=@(Need $mat 'optical_constants');if($optical.Count -ne 2){throw 'Exactly two frozen optical rows required:10.5/9um'}
  foreach($row in $rows){$matches=@($optical|Where-Object {$_.wavelength_um -eq $row.lambda_um});if($matches.Count -ne 1){throw 'No unique wavelength-bound material optical-constant row'};foreach($n in @('n_film','k_film','n_ti','k_ti')){Assert-EqualNumber $matches[0] $n (Num $row $n)}}
}
function Resources {
  $os=Get-CimInstance Win32_OperatingSystem;$cpu=(Get-CimInstance Win32_Processor|Measure-Object LoadPercentage -Average).Average
  if($null -eq $os.FreePhysicalMemory -or $null -eq $cpu){throw 'Resource counters unavailable'}
  $volume=@(Get-Volume -DriveLetter D -ErrorAction Stop);if($volume.Count -ne 1 -or !(Number $volume[0].SizeRemaining)){throw 'D scratch volume counters unavailable'}
  $health=[string]$volume[0].HealthStatus;$operations=@($volume[0].OperationalStatus|ForEach-Object {[string]$_});$diskOK=$health -eq 'Healthy' -and $operations.Count -gt 0 -and @($operations|Where-Object {$_ -ne 'OK'}).Count -eq 0
  $jobs=@(Get-CimInstance Win32_Process -Filter "Name LIKE '%comsol%' OR Name LIKE '%fdtd%' OR Name='mpiexec.exe' OR Name='mpirun.exe' OR Name='hydra_pmi_proxy.exe'"|Select-Object ProcessId,ParentProcessId,Name,CreationDate)
  return [pscustomobject]@{free_gib=[double]$os.FreePhysicalMemory/1MB;cpu_percent=[double]$cpu;D_free_gib=[double]$volume[0].SizeRemaining/1GB;D_health=$health;D_operational=$operations;D_healthy=$diskOK;compute_process_count=$jobs.Count;compute_processes=$jobs;verified_utc=[DateTime]::UtcNow.ToString('o')}
}
function Guard([long]$Rss,[double]$HostFree,[int]$PreviousHard,[int]$PreviousLow){
  if($Rss -lt 0 -or !(Number $HostFree) -or $HostFree -lt 0 -or $PreviousHard -lt 0 -or $PreviousLow -lt 0){throw 'Invalid resource sample'}
  $h=if($Rss -gt $HardGiB*1GB){$PreviousHard+1}else{0};$l=if($HostFree -lt $RuntimeFreeGiB){$PreviousLow+1}else{0};$reason=$null
  if($h -ge $GuardSamplesRequired){$reason='OWNED_RSS_HARD_LIMIT'}elseif($l -ge $GuardSamplesRequired){$reason='HOST_FREE_MEMORY_RESERVE'}
  return [pscustomobject]@{hard=$h;low=$l;reason=$reason;soft=($Rss -gt $SoftGiB*1GB)}
}

function Assert-ConditionalStartBudget {return}
function Stage-NativeJobs {
 $all=@(Get-CimInstance Win32_Process);$ids=@{};foreach($x in $all){if(([string]$x.CommandLine).IndexOf($stage,[StringComparison]::OrdinalIgnoreCase) -ge 0){$ids[[int]$x.ProcessId]=$true}}
 do{$more=$false;foreach($x in $all){if($ids.ContainsKey([int]$x.ParentProcessId) -and !$ids.ContainsKey([int]$x.ProcessId)){$ids[[int]$x.ProcessId]=$true;$more=$true}}}while($more)
 return @($all|Where-Object {$ids.ContainsKey([int]$_.ProcessId) -and $_.Name -match '^(?:comsol[^.]*|javaw?|mpiexec(?:\.hydra)?|hydra_pmi_proxy)\.exe$'})
}
function Shared-Admission($r,[int]$StageNativeCount){return $r.free_gib -ge 100 -and $r.cpu_percent -lt 75 -and $r.D_healthy -and $r.D_free_gib -ge 100 -and $StageNativeCount -eq 0}
function Wait-Resources {
 do{Assert-Window;if(Test-Path -LiteralPath (Join-Path $stage 'STOP_NEW_ADMISSION')){throw 'STOP_NEW_ADMISSION present'};Assert-Manifest;$r=Resources;$own=@(Stage-NativeJobs);$ok=Shared-Admission $r $own.Count;if(!$ok){Save-Status 'WAITING_FOR_RESOURCE_ADMISSION' $r;Start-Sleep -Seconds 30}}while(!$ok)
 Assert-Window;return $r
}
function Capture-OwnedRootIdentity($Process,[string]$ExpectedCommandFragment){
  $Process.Refresh();if($Process.HasExited){throw 'Owned CMD exited before identity capture'}
  $r=@(Get-CimInstance Win32_Process -Filter ('ProcessId='+$Process.Id));if($r.Count -ne 1 -or $null -eq $r[0].CreationDate -or [int]$r[0].ParentProcessId -ne [int]$PID){throw 'Cannot establish newly started owned CMD identity'}
  if([string]::IsNullOrEmpty([string]$r[0].CommandLine) -or ([string]$r[0].CommandLine).IndexOf($ExpectedCommandFragment,[StringComparison]::OrdinalIgnoreCase) -lt 0){throw 'Owned command line mismatch'}
  $Process.Refresh();if($Process.HasExited){throw 'Owned CMD exited during capture'}
  $script:ownedRootPid=[int]$r[0].ProcessId;$script:ownedRootCreated=([DateTime]$r[0].CreationDate).ToUniversalTime();$script:ownedIds=@{};$script:ownedIds[$script:ownedRootPid]=$script:ownedRootCreated
}
function Owned-Tree([int]$RootId){
  $all=@(Get-CimInstance Win32_Process);$root=@($all|Where-Object {[int]$_.ProcessId -eq $RootId})
  if($root.Count -ne 1 -or $null -eq $root[0].CreationDate -or ([DateTime]$root[0].CreationDate).ToUniversalTime() -ne $script:ownedRootCreated){throw 'Owned CMD missing or PID reused'}
  $ids=[Collections.Generic.HashSet[int]]::new();[void]$ids.Add($RootId);$birth=@{}
  foreach($x in $all){if($null -ne $x.CreationDate){$birth[[int]$x.ProcessId]=([DateTime]$x.CreationDate).ToUniversalTime()}}
  do{$added=$false;foreach($x in $all){$p=[int]$x.ParentProcessId;$i=[int]$x.ProcessId;if($ids.Contains($p) -and !$ids.Contains($i) -and $birth.ContainsKey($p) -and $birth.ContainsKey($i) -and $birth[$i] -ge $birth[$p]){[void]$ids.Add($i);$added=$true}}}while($added)
  $tree=@($all|Where-Object {$ids.Contains([int]$_.ProcessId)});foreach($x in $tree){$script:ownedIds[[int]$x.ProcessId]=([DateTime]$x.CreationDate).ToUniversalTime()};return $tree
}
function Identity-Present([int]$ProcessId,[DateTime]$CreatedUtc){$r=@(Get-CimInstance Win32_Process -Filter ('ProcessId='+$ProcessId) -ErrorAction Stop);if($r.Count -eq 0){return $false};if($r.Count -ne 1 -or $null -eq $r[0].CreationDate){throw 'Identity recheck unavailable'};return ([DateTime]$r[0].CreationDate).ToUniversalTime() -eq $CreatedUtc}
function Ensure-HandleApi {
  if($null -ne ('CodexCoatedOwnedHandles' -as [type])){return}
  Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class CodexCoatedOwnedHandles {
 [DllImport("kernel32.dll",SetLastError=true)] public static extern IntPtr OpenProcess(uint access,bool inherit,uint pid);
 [DllImport("kernel32.dll",SetLastError=true)] [return:MarshalAs(UnmanagedType.Bool)] public static extern bool TerminateProcess(IntPtr handle,uint code);
 [DllImport("kernel32.dll",SetLastError=true)] [return:MarshalAs(UnmanagedType.Bool)] public static extern bool CloseHandle(IntPtr handle);
}
'@ -ErrorAction Stop
}
function Stop-OwnedIdentities {
  if($script:ownedIds.Count -eq 0){return};Ensure-HandleApi;$fail=[Collections.Generic.List[string]]::new()
  foreach($e in @($script:ownedIds.GetEnumerator()|Sort-Object Value -Descending)){
    $handle=[IntPtr]::Zero;$errorText=$null
    try{
      # Pin the kernel object, then recheck the exact original CIM creation time.
      $handle=[CodexCoatedOwnedHandles]::OpenProcess([uint32]0x00101001,$false,[uint32]$e.Key)
      if($handle -eq [IntPtr]::Zero){$errorText='OpenProcess failed for owned PID '+$e.Key}
      elseif(Identity-Present ([int]$e.Key) ([DateTime]$e.Value)){if(![CodexCoatedOwnedHandles]::TerminateProcess($handle,0)){$errorText='TerminateProcess(handle) failed for owned PID '+$e.Key}}
    }catch{$errorText=$_.Exception.Message}
    finally{if($handle -ne [IntPtr]::Zero -and ![CodexCoatedOwnedHandles]::CloseHandle($handle)){$fail.Add('CloseHandle failed for owned PID '+$e.Key)}}
    if($null -ne $errorText){try{if(Identity-Present ([int]$e.Key) ([DateTime]$e.Value)){$fail.Add($errorText)}}catch{$fail.Add($errorText+'; identity recheck failed: '+$_.Exception.Message)}}
  }
  if($fail.Count){throw ($fail.ToArray() -join '; ')}
}
function Log-Tail([string]$Path){
  if(!(Test-Path -LiteralPath $Path -PathType Leaf)){return ''}
  $f=[IO.File]::Open($Path,[IO.FileMode]::Open,[IO.FileAccess]::Read,[IO.FileShare]::ReadWrite)
  try{$count=[int][Math]::Min(16384L,$f.Length);[void]$f.Seek(-$count,[IO.SeekOrigin]::End);$bytes=New-Object byte[] $count;$got=$f.Read($bytes,0,$count);return [Text.Encoding]::UTF8.GetString($bytes,0,$got)}finally{$f.Dispose()}
}
function Marker([string]$Path,[string]$Prefix){
  if(!(Test-Path -LiteralPath $Path) -or (Get-Item -LiteralPath $Path).Length -gt 67108864){throw 'Missing or oversized native stdout; manual review required'}
  $reader=[IO.StreamReader]::new($Path,[Text.Encoding]::UTF8);$matches=[Collections.Generic.List[string]]::new()
  try{while($null -ne ($line=$reader.ReadLine())){if($line.StartsWith($Prefix)){$matches.Add($line.Substring($Prefix.Length))}}}finally{$reader.Dispose()}
  if($matches.Count -ne 1){$failure=Native-Failure $Path;if($null -ne $failure){throw (Native-FailureMessage $failure)};throw ('Expected exactly one '+$Prefix+'; got '+$matches.Count)};return $matches[0]|ConvertFrom-Json
}
function Native-Failure([string]$Path){
  if(!(Test-Path -LiteralPath $Path -PathType Leaf)){return $null};if((Get-Item -LiteralPath $Path).Length -gt 67108864){throw 'Oversized native stdout; failure evidence requires manual review'}
  $reader=[IO.StreamReader]::new($Path,[Text.Encoding]::UTF8);$first=$null;$count=0
  try{while($null -ne ($line=$reader.ReadLine())){if($line.StartsWith('CODEX_FAILURE_JSON=')){$entry=$line.Substring('CODEX_FAILURE_JSON='.Length)|ConvertFrom-Json;if($null -eq $first){$first=$entry};$count++}}}finally{$reader.Dispose()}
  if($null -ne $first){$first|Add-Member -NotePropertyName native_failure_marker_count -NotePropertyValue $count -Force};return $first
}
function Native-FailureMessage($Failure){return ('Native failure in '+[string](Need $Failure 'failed_phase')+': '+[string](Need $Failure 'exception_type')+': '+[string](Need $Failure 'exception_message')+'; original CODEX_FAILURE_JSON and complete stdout retained')}
function Run-OwnedCmd([string]$Directory,[string]$CmdFile){
  $p=Start-Process -FilePath $env:ComSpec -ArgumentList @('/d','/s','/c',('"'+$CmdFile+'"')) -WorkingDirectory $Directory -PassThru -WindowStyle Hidden
  Capture-OwnedRootIdentity $p $CmdFile;$script:active.process_id=$p.Id;$hard=0;$low=0
  while(!$p.HasExited){
    Start-Sleep -Seconds 5;$p.Refresh();if($p.HasExited){break}
    try{$tree=@(Owned-Tree $p.Id)}catch{$p.Refresh();if($p.HasExited){break};throw}
    $rss=0L;foreach($x in $tree){$rss+=[long]$x.WorkingSetSize};$script:peakRss=[Math]::Max($script:peakRss,$rss);$script:active.owned_rss_bytes=$rss;$script:active.owned_process_ids=@($tree|ForEach-Object {[int]$_.ProcessId})
    $tail=Log-Tail (Join-Path $Directory 'stdout.log');$progress=@($tail -split "`n"|Where-Object {$_.StartsWith('CODEX_PROGRESS_JSON=')})
    if($progress.Count){try{$q=$progress[-1].Substring('CODEX_PROGRESS_JSON='.Length)|ConvertFrom-Json;$script:active.native_progress=$q;Write-Json (Join-Path $Directory 'progress.json') $q}catch{}}
    Assert-Manifest;if($Mode -eq 'SolveApproved' -and (Hash (Admission-Path)) -ne $admissionHash){throw 'Admission changed during physical attempt'};$r=Resources;if(!$r.D_healthy -or $r.D_free_gib -lt 100){throw 'D scratch volume unhealthy or free space below100GiB'};$g=Guard $rss $r.free_gib $hard $low;$hard=$g.hard;$low=$g.low;$script:active.hard_samples=$hard;$script:active.host_low_samples=$low
    if($null -ne $g.reason){Save-Status ('RESOURCE_STOP_'+$g.reason) $r $g.reason;Stop-OwnedIdentities;throw $g.reason}
    if($g.soft){Save-Status 'RUNNING_SOFT_MEMORY_WARNING' $r}else{Save-Status 'RUNNING' $r}
  }
  $exit=Join-Path $Directory 'exitcode.txt';if(!(Test-Path -LiteralPath $exit)){throw 'Owned CMD finished without exitcode artifact'}
  $script:ownedRootPid=$null;$script:ownedRootCreated=$null;$script:ownedIds=@{}
  return [int]([IO.File]::ReadAllText($exit,[Text.Encoding]::ASCII).Trim())
}
function Assert-Geometry($Geometry){
  if((Need $Geometry 'case_id') -ne $CaseId -or !(True (Need $Geometry 'geometry_verified')) -or (Need $Geometry 'geometry') -ne (Case-GeometryMarker) -or (Need $Geometry 'surface') -ne 'TI_IBC_WITH_SIO2_FLOOR_FILL' -or !(True (Need $Geometry 'not_conformal_wall_or_land_coating')) -or (Need $Geometry 'proof_schema') -ne 'ANGLE_BAND_SOURCE_POWER_V2'){throw 'Fixed case-bound G99 floor-fill geometry/source proof not bound/verified'}
  foreach($x in @(@('wavelength_um','lambda_um'),@('theta_deg','theta_deg'),@('phi_global_deg','phi_deg'))){Assert-EqualNumber $Geometry $x[0] (Num $case $x[1])}
  $normal=Need $Geometry 'normal_incidence_only';if($normal -isnot [bool] -or $normal -ne ((Num $case 'theta_deg') -eq 0)){throw 'Normal/oblique incidence declaration is not case-bound'}
  $thickness=Num $case 'film_um';Assert-EqualNumber $Geometry 'floor_fill_thickness_um' $thickness;Assert-EqualNumber $Geometry 'bottom_flat_to_flat_um' (Num $case 'bottom_flat_to_flat_um')
  Assert-ProfileBinding $Geometry $case
  Assert-EqualNumber $Geometry 'opening_flat_to_flat_um' (Num $case 'opening_flat_to_flat_um')
  $measures=Analytic-Geometry $case
  foreach($n in $measures.Keys){if([Math]::Abs((Num $Geometry $n)/$measures[$n]-1) -gt 1e-6){throw ('Measured floor-fill geometry differs from independent analytic dimensions: '+$n)}}
  $expect=Need $case 'expected_geometry';$mapping=@{film_top_apothem_um=@('film_top_apothem_um',1.);film_volume_m3=@('film_volume_um3',1e-18);air_volume_m3=@('air_volume_um3',1e-18);total_computational_volume_m3=@('total_dielectric_volume_um3',1e-18);film_air_interface_area_m2=@('film_air_interface_area_um2',1e-12);top_area_m2=@('top_area_um2',1e-12)}
  foreach($n in $mapping.Keys){$bound=$mapping[$n];if([Math]::Abs((Num $Geometry $n)/((Num $expect $bound[0])*[double]$bound[1])-1) -gt 1e-6){throw ('Measured geometry differs from frozen per-case expected geometry: '+$n)}}
  $fd=Num $Geometry 'film_domain';$ad=Num $Geometry 'air_domain';if($fd -le 0 -or $ad -le 0 -or $fd -eq $ad -or $fd -ne [Math]::Floor($fd) -or $ad -ne [Math]::Floor($ad)){throw 'Distinct actual film/air domains required'}
  if(@(Need $Geometry 'periodic_side_ids').Count -ne 6 -or @(Need $Geometry 'film_interface_faces').Count -ne 1 -or @(Need $Geometry 'Ti_film_side_faces').Count -ne 7 -or @(Need $Geometry 'Ti_air_side_faces').Count -ne (6*(Case-Contract $CaseId).profile_segments+1)){throw 'Fixed film interface/Ti adjacency inventory mismatch'}
  $inc=Incidence $case;Assert-Vector (Need $Geometry 'global_E_direction') $inc.E_in 1e-12 'incoming case-bound true p/s polarization direction'
  foreach($n in @('height_um','air_cap_um','top_plane_z_um')){$wanted=if($n -eq 'height_um'){Num $case 'height_um'}elseif($n -eq 'air_cap_um'){10.5}else{(Num $case 'height_um')+10.5};Assert-EqualNumber $Geometry $n $wanted};Assert-Vector (Need $Geometry 'port_reference_xyz_m') @(40e-6,(40e-6/[Math]::Sqrt(3)),(((Num $case 'height_um')+10.5)*1e-6)) 1e-11 'fixed top-port reference'
}

function Assert-MeshProfile($Mesh,$Row){
 $short=(Need $Row 'mesh') -eq 'm4_air29_wft3_film02_short9';$variant=if($short){'short9_film02_Ti075_edge0225'}else{'original_m4_air29_wft3'}
 if((Need $Mesh 'mesh') -ne (Need $Row 'mesh') -or (Need $Mesh 'mesh_local_variant') -ne $variant -or !(True (Need $Mesh 'mesh_local_top_boundary_size_readback_verified')) -or $null -ne (Need $Mesh 'mesh_hLocalTop_um') -or !(True (Need $Mesh 'case_bound_mesh_size_readback_verified'))){throw 'Case-bound dynamic exact mesh-size readback missing'}
 $e=Need $Row 'expected_geometry';foreach($x in @(@('mesh_hAir_um','air_hmax_um'),@('mesh_hTi_um','Ti_hmax_um'),@('mesh_hEdge_um','edge_hmax_um'),@('mesh_hFilm_um','film_hmax_um'),@('film_hmax_um','film_hmax_um'))){Assert-EqualNumber $Mesh $x[0] (Num $e $x[1])}
 $p=Need $Mesh 'loss_region_mesh_admission';if((Need $p 'profile') -ne $Row.mesh -or !(True (Need $p 'actual_profile_size_readback_verified')) -or !(True (Need $p 'conditional_loss_mesh_admission')) -or (Need $p 'conditional9') -isnot [bool] -or [bool]$p.conditional9 -ne $short){throw 'Actual profile/loss-region mesh admission missing'}
 foreach($x in @(@('actual_air_hmax_um','air_hmax_um'),@('actual_Ti_boundary_hmax_um','Ti_hmax_um'),@('actual_edge_hmax_um','edge_hmax_um'),@('actual_film_hmax_um','film_hmax_um'))){Assert-EqualNumber $p $x[0] (Num $e $x[1])}
 $decay=(Num $Row 'lambda_um')/(4*[Math]::PI*(Num $Row 'k_film'));Assert-EqualNumber $p 'film_intensity_decay_length_um' $decay;Assert-EqualNumber $p 'film_hmax_over_intensity_decay' ((Num $e 'film_hmax_um')/$decay);Assert-EqualNumber $p 'conditional9_maximum_film_hmax_over_intensity_decay' .5
 if($short -and ((Num $Row 'lambda_um') -ne 9 -or (Num $p 'film_hmax_over_intensity_decay') -gt .5)){throw 'Conditional9 film attenuation-length mesh guard failed'}
}
function Assert-Mesh($Mesh,$Geometry){
  Assert-Geometry $Geometry;Assert-OOC (Need $Mesh 'linear_solver_override')
  Assert-PresolvePortProof (Need $Mesh 'presolve_port_wavevector_proof')
  if((Need $Mesh 'case_id') -ne $CaseId -or (Need $Mesh 'state') -ne 'MESH_PREFLIGHT_ONLY' -or (Need $Mesh 'solver_executed') -isnot [bool] -or $Mesh.solver_executed -or (Need $Mesh 'mesh') -ne $case.mesh){throw 'Actual case-bound m4_air29_wft3 mesh-only receipt required'}
  foreach($n in @('geometry_verified','periodic_mesh_verified','second_order_curl_elements_requested')){if(!(True (Need $Mesh $n))){throw ('Mesh verification failed: '+$n)}}
  if((Need $Mesh 'reported_dof_kind') -ne 'actual' -or (Need $Mesh 'electric_curl_element_order') -ne '2'){throw 'Actual DOF and quadratic curl elements required'}
  if((Need $Mesh 'mesh_fingerprint_sha256') -notmatch '^[a-f0-9]{64}$' -or (Need $Mesh 'mesh_fingerprint_scheme') -ne 'SHA256_RAW_COORDINATE_BITS_TET_TRI_CONNECTIVITY_ENTITY_IDS_V1'){throw 'Missing exact-coordinate/connectivity actual mesh fingerprint'}
  $resourceBounds=Need $script:plan 'resource_plan'
  $actualDof=Num $Mesh 'degrees_of_freedom';$actualTet=Num $Mesh 'tetrahedra';$actualQuality=Num $Mesh 'min_tetrahedron_quality'
  foreach($x in @(@('degrees_of_freedom',$actualDof),@('tetrahedra',$actualTet))){$v=[double]$x[1];if($v -le 0 -or $v -ne [Math]::Floor($v)){throw ('Actual mesh admission inventory is not a positive integer: '+$x[0])}}
  if($actualQuality -lt 0 -or $actualQuality -gt 1){throw 'Actual tetrahedron quality is outside the finite physical range0..1'}
  $modes=Need $case 'expected_modes';if((Num $Mesh 'propagating_spatial_orders') -ne (Num $modes 'spatial_propagating_orders') -or (Num $Mesh 'polarization_channels') -ne (Num $modes 'propagating_polarization_channels') -or (Num $Mesh 'total_port_channels_including_evanescent') -ne (Num $modes 'all_polarization_channels')){throw 'Case-bound incidence-shifted propagating/padded channel inventory mismatch'}
  if((Need $Mesh 'proof_schema') -ne 'ANGLE_BAND_SOURCE_POWER_V2' -or (Need $Mesh 'polarization') -ne $case.polarization){throw 'Mesh source proof schema/polarization mismatch'}
  foreach($x in @(@('wavelength_um','lambda_um'),@('theta_deg','theta_deg'),@('phi_global_deg','phi_deg'))){Assert-EqualNumber $Mesh $x[0] (Num $case $x[1])}
  foreach($n in @('requested_incidence_verified','global_phi_verified','side_floquet_source_verified')){if(!(True (Need $Mesh $n))){throw ('Actual oblique-source verification failed before admission: '+$n)}}
  $inc=Incidence $case;Assert-Vector (Need $Mesh 'expected_requested_k_in_rad_per_m') $inc.k_in .1 'mesh requested physical k-in';Assert-Vector (Need $Mesh 'actual_k_in_rad_per_m') $inc.k_in .1 'mesh actual physical k-in'
  Assert-EqualNumber $Mesh 'floor_fill_thickness_um' (Num $case 'film_um');Assert-EqualNumber $Mesh 'film_hmax_um' (Num (Need $case 'expected_geometry') 'film_hmax_um')
  foreach($x in @(@('mesh_hAir_um','air_hmax_um'),@('mesh_hTi_um','Ti_hmax_um'),@('mesh_hEdge_um','edge_hmax_um'),@('mesh_hFilm_um','film_hmax_um'))){Assert-EqualNumber $Mesh $x[0] (Num (Need $case 'expected_geometry') $x[1])}
  if([Math]::Abs((Num $Mesh 'film_volume_m3')/(Num $Geometry 'film_volume_m3')-1) -gt 1e-6 -or $Mesh.film_domain -ne $Geometry.film_domain -or $Mesh.air_domain -ne $Geometry.air_domain){throw 'Actual mesh material partition mismatch'}
  $pairs=@(Need $Mesh 'periodic_mesh_pairs');if($pairs.Count -ne 3){throw 'Three matched homogeneous-air periodic mesh pairs required'}
  foreach($p in $pairs){if(!(True (Need $p 'canonical_triangle_multiset_equal')) -or (Num $p 'source_triangles') -le 0 -or (Num $p 'source_triangles') -ne (Num $p 'destination_triangles')){throw 'Actual periodic triangle multisets mismatch'}}
  # Only after every geometry/source/periodic/OOC/fingerprint invariant above is
  # verified may a bounded mesh-only resource/quality rejection be typed. No
  # missing, nonfinite, nonpositive or noninteger evidence can be classified as
  # a skippable resource issue, and SolveApproved never emits a skip allowance.
  Assert-MeshProfile $Mesh $case
  $kind=$null;$reason=$null
  if($actualDof -gt (Num $resourceBounds 'actual_DOF_max')){$kind='DOF_CAP';$reason='Actual mesh admission inventory exceeds bounds: degrees_of_freedom'}
  elseif($actualTet -gt (Num $resourceBounds 'tetrahedra_max')){$kind='TET_CAP';$reason='Actual mesh admission inventory exceeds bounds: tetrahedra'}
  elseif($actualQuality -lt (Num $resourceBounds 'min_tetrahedron_quality')){$kind='MIN_QUALITY';$reason='Actual tetrahedron quality below0.2'}
  if($null -ne $kind){
    if($Mode -eq 'MeshOnly'){$script:meshPrecheckRejection=[ordered]@{schema=1;kind=$kind;case_id=$CaseId;mode=$Mode;degrees_of_freedom=$actualDof;actual_DOF_max=Num $resourceBounds 'actual_DOF_max';tetrahedra=$actualTet;actual_tetrahedra_max=2500000;min_tetrahedron_quality=$actualQuality;min_tetrahedron_quality_required=Num $resourceBounds 'min_tetrahedron_quality';geometry_source_periodic_OOC_and_fingerprint_checks_passed=$true;solver_executed=$false;automatic_coarsening_allowed=$false;reason=$reason}}
    throw $reason
  }
}
function Assert-Receipt([string]$Directory){
  $r=Json (Join-Path $Directory 'case_receipt.json')
  foreach($b in @(@('core_template_sha256','CavityFloorFilm.java'),@('controller_sha256','run_floorfilm_controller.ps1'),@('plan_sha256','plan.json'),@('material_sha256','material.json'))){if((Need $r $b[0]) -ne (Hash (Join-Path $stage $b[1]))){throw 'Receipt source hash no longer matches'}}
  if($r.case_id -ne $CaseId -or $r.manifest_sha256 -ne $manifestHash -or $r.generated_source_sha256 -ne (Hash (Join-Path $Directory 'CaseEntry.java')) -or $r.compiled_case_sha256 -ne (Hash (Join-Path $Directory 'CaseEntry.class'))){throw 'Case source/class receipt binding failed'}
  if((Need $r 'queue_window_sha256') -ne $windowHash -or (Need $r 'linear_solver_override_sha256') -ne (Hash (Join-Path $Directory 'linear_solver_override.json'))){throw 'Receipt queue-window/OOC readback binding changed'}
  $proofPath=Join-Path $Directory 'presolve_port_wavevectors.json';if((Need $r 'presolve_port_wavevector_proof_sha256') -ne (Hash $proofPath)){throw 'Receipt presolve actual-port proof hash changed'};$proof=Json $proofPath;Assert-PresolvePortProof $proof
  $mesh=Json (Join-Path $Directory 'mesh.json');if(($proof|ConvertTo-Json -Depth 30 -Compress) -ne ((Need $mesh 'presolve_port_wavevector_proof')|ConvertTo-Json -Depth 30 -Compress)){throw 'Independent presolve marker proof differs from saved mesh inventory'}
  Assert-OOC (Json (Join-Path $Directory 'linear_solver_override.json'));return $r
}
function Assert-Admission {
  Assert-PairAdmission
  $path=Admission-Path;$h=Hash $path;$a=Json $path
  if($null -ne $script:admissionHash -and $h -ne $script:admissionHash){throw 'Admission changed during this attempt'}
  if($a.schema -ne 1 -or $a.purpose -ne 'FLOORFILM_CAVITY_MESH_ADMISSION' -or $a.stage_root -ne $stage -or [string]::IsNullOrWhiteSpace([string]$a.reviewed_by) -or [string]::IsNullOrWhiteSpace([string]$a.reviewed_utc)){throw 'Missing explicit stage-bound mesh review provenance'}
  if((Need $a 'queue_window_sha256') -ne $windowHash){throw 'Admission belongs to another queue window'}
  $rows=@($a.entries|Where-Object {$_.case_id -eq $CaseId});if(@($a.entries).Count -ne 1 -or $rows.Count -ne 1 -or !(True (Need $rows[0] 'allow_solve'))){throw 'No unique case-isolated reviewed solve admission'};$e=$rows[0]
  foreach($b in @(@('core_sha256','CavityFloorFilm.java'),@('plan_sha256','plan.json'),@('controller_sha256','run_floorfilm_controller.ps1'))){if((Need $e $b[0]) -ne (Hash (Join-Path $stage $b[1]))){throw 'Admission source hash mismatch'}}
  if((Need $e 'manifest_sha256') -ne $manifestHash){throw 'Admission manifest hash mismatch'};if((Need $e 'fresh_flat_admission_path') -ne (Flat-Path $CaseId) -or (Need $e 'fresh_flat_admission_sha256') -ne (Assert-ConditionalSourceAdmission)){throw 'Solve admission not bound to immutable fresh native source control'}
  $rel=[string](Need $e 'mesh_relative_directory');$dir=[IO.Path]::GetFullPath((Join-Path $stage $rel));$base=(Join-Path $stage $CaseId)+'\'
  if([IO.Path]::IsPathRooted($rel) -or !$dir.StartsWith($base,[StringComparison]::OrdinalIgnoreCase) -or [IO.Path]::GetFileName($dir) -notmatch '^mesh_preflight(?:_[A-Za-z0-9_-]+)?$'){throw 'Admission path is not a case-owned mesh preflight'}
  foreach($b in @(@('mesh_sha256','mesh.json'),@('geometry_sha256','geometry.json'),@('receipt_sha256','case_receipt.json'),@('mesh_model_sha256','mesh_only.mph'),@('presolve_port_wavevector_proof_sha256','presolve_port_wavevectors.json'))){if((Need $e $b[0]) -ne (Hash (Join-Path $dir $b[1]))){throw 'Reviewed mesh artifacts changed'}}
  [void](Assert-Receipt $dir);$mesh=Json (Join-Path $dir 'mesh.json');Assert-Mesh $mesh (Json (Join-Path $dir 'geometry.json'))
  $script:admissionHash=$h;return $mesh
}

function Assert-ConditionalSourceAdmission {
 $path=Flat-Path $CaseId;$a=Json $path
 if((Need $a 'state') -ne (Flat-State $CaseId) -or (Need $a 'phase') -ne (Source-Phase $CaseId) -or (Need $a 'manifest_sha256') -ne $manifestHash -or (Need $a 'queue_window_sha256') -ne $windowHash -or @(Need $a 'points').Count -ne (Flat-Count $CaseId)){throw 'Immutable phase-bound fresh native source admission missing'}
 $points=@($a.points|Where-Object {$_.approved_case_id -eq $CaseId});if($points.Count -ne 1){throw 'No unique exact-case fresh source control'};$p=$points[0]
 if($p.state -ne 'PASS' -or $p.polarization -ne $case.polarization -or $p.wavelength_um -ne $case.lambda_um -or $p.theta_deg -ne $case.theta_deg -or $p.phi_global_deg -ne $case.phi_deg){throw 'Actual phase-flat parameters differ from requested case'}
 foreach($h in $p.artifacts_sha256.PSObject.Properties){if([IO.Path]::GetFileName($h.Name) -ne $h.Name -or (Hash (Join-Path $p.directory $h.Name)) -ne $h.Value){throw 'Fresh phase-flat source artifact changed'}};return (Hash $path)
}
function Java-Double([double]$Value){return $Value.ToString('R',[Globalization.CultureInfo]::InvariantCulture)}
function Write-CaseEntry([string]$Directory,[bool]$MeshOnly,$ReviewedMesh){
  $core=Join-Path $stage 'CavityFloorFilm.java';$text=[IO.File]::ReadAllText($core,[Text.Encoding]::UTF8);$header='public class CavityFloorFilm {'
  if([regex]::Matches($text,[regex]::Escape($header)).Count -ne 1){throw 'Expected exactly one CavityFloorFilm class declaration'}
  $literal=$Directory.Replace('\','\\');$flag=if($MeshOnly){'true'}else{'false'};$review=''
  if(!$MeshOnly){$review='reviewedDof='+([long]$ReviewedMesh.degrees_of_freedom)+'L; reviewedTet='+([long]$ReviewedMesh.tetrahedra)+'L; reviewedMeshSha="'+$ReviewedMesh.mesh_fingerprint_sha256+'"; '}
  $javaArgs=@('n_film','k_film')|ForEach-Object {Java-Double ([double]$case.$_)}
  $insert='public class CaseEntry { public static void main(String[] args) throws Exception { '+$review+'run("'+$literal+'","'+$CaseId+'",'+($javaArgs -join ',')+',"'+$case.mesh+'",'+$flag+'); }'
  [IO.File]::WriteAllText((Join-Path $Directory 'CaseEntry.java'),$text.Replace($header,$insert),[Text.UTF8Encoding]::new($false))
  $receipt=[ordered]@{schema=1;case_id=$CaseId;mode=$Mode;attempt_id=$AttemptId;mesh_only=$MeshOnly;parameters=$case;manifest_sha256=$manifestHash;admission_sha256=$admissionHash;queue_window_sha256=$windowHash;linear_solver_override_sha256=$null;core_template_sha256=Hash $core;controller_sha256=Hash (Join-Path $stage 'run_floorfilm_controller.ps1');plan_sha256=Hash (Join-Path $stage 'plan.json');material_sha256=Hash (Join-Path $stage 'material.json');generated_source_sha256=Hash (Join-Path $Directory 'CaseEntry.java');compiled_case_sha256=$null;reviewed_actual_DOF=$(if($MeshOnly){$null}else{$ReviewedMesh.degrees_of_freedom});reviewed_actual_tetrahedra=$(if($MeshOnly){$null}else{$ReviewedMesh.tetrahedra});reviewed_actual_mesh_fingerprint_sha256=$(if($MeshOnly){$null}else{$ReviewedMesh.mesh_fingerprint_sha256})}
  New-Json (Join-Path $Directory 'case_receipt.json') $receipt
  $scratch=Join-Path (Join-Path (Join-Path $stage 'scratch') $CaseId) (($Mode.ToLowerInvariant())+'_'+$AttemptId)
  if(Test-Path -LiteralPath $scratch){throw 'Scratch already exists; no overwrite/retry'};New-Item -ItemType Directory -Path $scratch -ErrorAction Stop|Out-Null
  $hashCmd="`$ErrorActionPreference='Stop';[IO.File]::WriteAllText('$Directory\compiled_case.sha256',(Get-FileHash -LiteralPath '$Directory\CaseEntry.class' -Algorithm SHA256).Hash.ToLowerInvariant(),[Text.Encoding]::ASCII)"
  $encoded=[Convert]::ToBase64String([Text.Encoding]::Unicode.GetBytes($hashCmd));$psExe=Join-Path $env:WINDIR 'System32\WindowsPowerShell\v1.0\powershell.exe'
  $compiler=Join-Path $ComsolBin 'comsolcompile.exe';$batch=Join-Path $ComsolBin 'comsolbatch.exe'
  $cmd="@echo off`r`nset TEMP=$scratch`r`nset TMP=$scratch`r`ncd /d `"$Directory`"`r`n`"$compiler`" `"$Directory\CaseEntry.java`" > `"$Directory\compile.stdout.log`" 2>&1`r`nset CODEX_COMPILE_EXIT=%ERRORLEVEL%`r`n> `"$Directory\compile_exitcode.txt`" echo %CODEX_COMPILE_EXIT%`r`nif not %CODEX_COMPILE_EXIT%==0 goto compile_failed`r`n"
  $cmd+="if not exist `"$Directory\CaseEntry.class`" goto compile_failed`r`nfindstr /C:`"Compilation failed`" `"$Directory\compile.stdout.log`" > nul`r`nif not errorlevel 1 goto compile_failed`r`n"
  $cmd+="`"$psExe`" -NoProfile -NonInteractive -EncodedCommand $encoded`r`nif not exist `"$Directory\compiled_case.sha256`" goto compile_failed`r`n"
  $cmd+="`"$batch`" -np $Cores -tmpdir `"$scratch`" -inputfile `"$Directory\CaseEntry.class`" -batchlog `"$Directory\solver.log`" -autosave off > `"$Directory\stdout.log`" 2>&1`r`n> `"$Directory\exitcode.txt`" echo %ERRORLEVEL%`r`nexit /b`r`n:compile_failed`r`n> `"$Directory\exitcode.txt`" echo 2`r`nexit /b 2`r`n"
  [IO.File]::WriteAllText((Join-Path $Directory 'run_case.cmd'),$cmd,[Text.Encoding]::ASCII);return $receipt
}
function Exploration-Disposition($Checks){
  $allowed=@('native_R_orders_vs_net_flux','native_flux_vs_total_loss_absolute','native_flux_vs_total_loss_relative','mode_vs_net','net_vs_loss_abs','net_vs_loss_relative')
  $failed=@($Checks.Keys|Where-Object {!(True $Checks[$_])});$warnings=@($failed|Where-Object {$_ -in $allowed});$hard=@($failed|Where-Object {$_ -notin $allowed})
  return [ordered]@{accepted=($failed.Count -eq 0);failed_checks=$failed;warning_checks=$warnings;hard_failed_checks=$hard;exploration_continue_allowed=($hard.Count -eq 0);disposition=$(if($hard.Count){'STOP_HARD_FAILURE'}elseif($warnings.Count){'COMPLETE_WITH_FLUX_WARNING'}else{'COMPLETE_EXPLORATORY_POINT_ACCEPTED'});formal_qualification=$false;mesh_convergence_established=$false;warning_policy='ONLY_SIX_ENUMERATED_NET_FLUX_CHECKS_ARE_NONBLOCKING_ORIGINAL_FLAGS_AND_THRESHOLDS_UNCHANGED'}
}
function Evaluate-Metrics($m,$ActualMesh,$ActualGeometry){
  if($m.case_id -ne $CaseId -or !(True (Need $m 'solver_executed')) -or !(True (Need $m 'solver_completed'))){throw 'Metrics require the exact completed physical case'}
  Assert-PresolvePortProof (Need $m 'presolve_port_wavevector_proof')
  Assert-MeshProfile $ActualMesh $case
  foreach($n in @('case_bound_mesh_size_readback_verified')){if(!(True (Need $m $n))){throw 'Solved mesh-size readback missing'}}
  if(((Need $m 'loss_region_mesh_admission')|ConvertTo-Json -Depth 20 -Compress) -ne ((Need $ActualMesh 'loss_region_mesh_admission')|ConvertTo-Json -Depth 20 -Compress)){throw 'Solved loss-mesh proof differs from admitted actual mesh'}
  if(((Need $m 'presolve_port_wavevector_proof')|ConvertTo-Json -Depth 30 -Compress) -ne ((Need $ActualMesh 'presolve_port_wavevector_proof')|ConvertTo-Json -Depth 30 -Compress)){throw 'Physical metrics do not preserve the actual pre-run complete-port wavevector proof'}
  foreach($x in @(@('wavelength_um','lambda_um'),@('theta_deg','theta_deg'),@('phi_global_deg','phi_deg'),@('nti','n_ti'),@('kti','k_ti'))){Assert-EqualNumber $m $x[0] (Num $case $x[1])}
  foreach($x in @(@('T',0),@('Pin_specified_W',1),@('film_hmax_um',(Num (Need $case 'expected_geometry') 'film_hmax_um')))){Assert-EqualNumber $m $x[0] $x[1]}
  if((Need $m 'proof_schema') -ne 'ANGLE_BAND_SOURCE_POWER_V2'){throw 'Metrics must preserve the actual oblique source/power-proof schema'}
  Assert-EqualNumber $m 'bottom_flat_to_flat_um' (Num $case 'bottom_flat_to_flat_um');Assert-EqualNumber $m 'floor_fill_thickness_um' (Num $case 'film_um');Assert-EqualNumber $m 'height_um' (Num $case 'height_um');Assert-EqualNumber $m 'air_cap_um' 10.5;Assert-EqualNumber $m 'top_plane_z_um' ((Num $case 'height_um')+10.5)
  if($m.polarization -ne $case.polarization -or $m.mesh -ne $case.mesh -or $m.geometry -ne (Case-GeometryMarker) -or $m.surface -ne 'TI_IBC_WITH_SIO2_FLOOR_FILL' -or $m.R_A_power_normalization -ne 'specified_1W_unchanged_not_fitted_Pin' -or !(True (Need $m 'not_conformal_wall_or_land_coating')) -or !(True (Need $m 'planar_TMM_is_not_cavity_acceptance_reference'))){throw 'Fixed case-bound thickness/polarization or unchanged1W normalization not preserved'}
  Assert-EqualNumber $m 'nfilm' (Num $case 'n_film');Assert-EqualNumber $m 'kfilm' (Num $case 'k_film')
  foreach($x in @(@('mesh_hAir_um',2.9),@('mesh_hTi_um',.875),@('mesh_hEdge_um',.2625),@('mesh_hFilm_um',.5))){Assert-EqualNumber $m $x[0] $x[1]}
  if((Need $m 'mesh_fingerprint_sha256') -ne $ActualMesh.mesh_fingerprint_sha256 -or (Need $m 'reviewed_mesh_fingerprint_sha256') -ne $ActualMesh.mesh_fingerprint_sha256 -or (Num $m 'degrees_of_freedom') -ne $ActualMesh.degrees_of_freedom -or (Num $m 'tetrahedra') -ne $ActualMesh.tetrahedra -or [Math]::Abs((Num $m 'film_volume_m3')/(Num $ActualGeometry 'film_volume_m3')-1) -gt 1e-6){throw 'Metrics not bound to actual reviewed mesh and material volume'}
  Assert-OOC (Need $m 'linear_solver_override');if(!(True (Need $m 'linear_solver_override_reverified_before_physics'))){throw 'Production OOC not strongly reverified immediately before runAll'}
  foreach($n in @('requested_incidence_verified','global_phi_verified','side_floquet_source_verified')){if(!(True (Need $m $n))){throw ('Physical angular-source proof failed: '+$n)}}
  $expectedModes=Need $case 'expected_modes';foreach($x in @(@('propagating_spatial_orders','spatial_propagating_orders'),@('polarization_channels','propagating_polarization_channels'),@('total_port_channels_including_evanescent','all_polarization_channels'))){Assert-EqualNumber $m $x[0] (Num $expectedModes $x[1])}
  $checks=Need $m 'checks'
  $required=@('geometry_verified','periodic_mesh_verified','all_propagating_channels_present','actual_wavevectors_match_independent_lattice','incident_wavevector_points_down','finite','incident_power','R_orders_vs_net_flux','periodic_side_power','flux_vs_total_loss_absolute','flux_vs_total_loss_relative','passive','explicit_global_source_binding','normalized_source_mode_global_direction','internal_E_fit_finite','internal_E_fit_condition','internal_E_fit_relative_L2_residual','internal_E_fit_all_LOO_Pin','internal_E_fit_cut_areas','film_parameter_binding','explicit_film_passive','film_manual_vs_builtin_Qh_absolute','film_manual_vs_builtin_Qh_relative','Ti_H_vs_E_loss_absolute','Ti_H_vs_E_loss_relative','port_vs_total_loss_absolute','port_vs_total_loss_relative','loss_partition_complete','total_loss_passive','requested_global_incidence','side_floquet_source_binding')
  if(@($checks.PSObject.Properties).Count -ne $required.Count){throw 'Native boolean check inventory differs from the frozen floor-fill contract'}
  foreach($n in $required){if((Need $checks $n) -isnot [bool]){throw ('Native check must be an explicit boolean: '+$n)}}
  $test=[ordered]@{};foreach($p in $checks.PSObject.Properties){$test['native_'+$p.Name]=True $p.Value}
  $nf=Num $case 'n_film';$kf=Num $case 'k_film';$test.material_permittivity=([Math]::Abs((Num $m 'eps_f_real')-($nf*$nf-$kf*$kf)) -le 1e-10 -and [Math]::Abs((Num $m 'eps_f_imag')+2*$nf*$kf) -le 1e-10 -and (True (Need $m 'film_parameter_binding_verified')))
  $r=Num $m 'R_total';$af=Num $m 'A_flux';$film=Num $m 'A_film';$ti=Num $m 'A_Ti_H';$tiE=Num $m 'A_Ti_E';$loss=Num $m 'A_loss_total';$pin=Num $m 'Pin_independent_internal_E_fit_W';$top=Num $m 'Pout_top_net_W'
  $port=1-$r;$gap=[Math]::Abs($port-$loss);$netgap=[Math]::Abs($af-$loss);$relative=$gap/[Math]::Max([Math]::Abs($port),1e-12);$netRelative=$netgap/[Math]::Max([Math]::Abs($af),1e-12)
  $tiGap=[Math]::Abs($tiE-$ti);$tiRelative=$tiGap/[Math]::Max([Math]::Abs($ti),1e-12);$filmQh=Num $m 'A_film_builtin_Qh';$filmGap=[Math]::Abs($film-$filmQh);$filmRelative=$filmGap/[Math]::Max([Math]::Abs($film),1e-12)
  $modeGap=[Math]::Abs($r-(1+$top))
  $test.loss_partition=([Math]::Abs($loss-$film-$ti) -le 1e-10);$test.port_A=([Math]::Abs((Num $m 'A_port')-$port) -le 1e-10);$test.top_flux_definition=([Math]::Abs($af+$top) -le 1e-10 -and [Math]::Abs((Num $m 'R_total_net_flux')-1-$top) -le 1e-10)
  $test.Ti_surface_partition=([Math]::Abs((Num $m 'A_floor')+(Num $m 'A_walls')+(Num $m 'A_land')-$ti) -le 1e-8 -and [Math]::Abs((Num $m 'A_Ti_film_side')+(Num $m 'A_Ti_air_side')-$ti) -le 1e-8)
  $test.mode_vs_net=($modeGap -le .005);$test.port_vs_loss_abs=($gap -le .005);$test.port_vs_loss_relative=($relative -le .02);$test.net_vs_loss_abs=($netgap -le .005);$test.net_vs_loss_relative=($netRelative -le .02)
  $test.Ti_E_vs_H_abs=($tiGap -le .005);$test.Ti_E_vs_H_relative=($tiRelative -le .02);$test.film_manual_vs_Qh_abs=($filmGap -le .005);$test.film_manual_vs_Qh_relative=($filmRelative -le .02)
  $test.incident_power=([Math]::Abs($pin-1) -le .005);$test.incident_error_record=([Math]::Abs((Num $m 'incident_power_relative_error')-[Math]::Abs($pin-1)) -le 1e-10)
  $test.passivity=($r -ge -.002 -and $r -le 1.002 -and $port -ge -.002 -and $port -le 1.002 -and $film -ge -1e-10 -and $ti -ge -1e-10 -and $loss -ge -.002 -and $loss -le 1.002 -and $af -ge -.002 -and $af -le 1.002)
  $fit=Need $m 'internal_E_power_fit';$incidence=Incidence $case;$k0=$incidence.k0;$ct=Assert-FitDirection $m $fit $case
  $test.source_fit=((Need $fit 'method') -eq 'FIXED_KZ_BLOCH_DEPHASED_FIVE_PLANE_TANGENTIAL_SP_COMPLEX_LS_V2' -and (Num $fit 'maximum_full_or_LOO_design_condition_2') -le 10 -and (Num $fit 'relative_L2_E_fit_residual') -le .005 -and (Num $fit 'LOO_max_Pin_relative_error_vs_specified_1W') -le .005 -and (True (Need $fit 'all_cut_areas_match_geometry')) -and [Math]::Abs((Num $fit 'Pin_W')-$pin) -le 1e-10 -and [Math]::Abs((Num $fit 'kz_fixed_rad_per_m')-$k0*$ct) -le .1 -and [Math]::Abs((Num $fit 'cos_theta')-$ct) -le 1e-12 -and [Math]::Abs((Num $fit 'power_admittance_p_over_vacuum')-1/$ct) -le 1e-12 -and [Math]::Abs((Num $fit 'power_admittance_s_over_vacuum')-$ct) -le 1e-12)
  Assert-Vector (Need $fit 'actual_k_in_rad_per_m') $incidence.k_in .1 'Bloch-dephased fit actual incoming wavevector'
  $loo=@(Need $fit 'leave_one_out');$test.source_LOO=($loo.Count -eq 5);foreach($x in $loo){$test.source_LOO=$test.source_LOO -and [Math]::Abs((Num $x 'Pin_W')-1) -le .005 -and (Num $x 'design_condition_2') -le 10}
  $planes=@(Need $fit 'plane_data');$fractions=@(.1,.25,.5,.75,.9);$area=2*[Math]::Sqrt(3)*(40e-6)*(40e-6);$test.five_fixed_source_planes=($planes.Count -eq 5)
  if($planes.Count -eq 5){for($j=0;$j -lt 5;$j++){$z=(Num $case 'height_um')*1e-6+$fractions[$j]*10.5e-6;$p=$planes[$j];$test.five_fixed_source_planes=$test.five_fixed_source_planes -and [Math]::Abs((Num $p 'cap_fraction')-$fractions[$j]) -le 1e-12 -and [Math]::Abs((Num $p 'z_m')-$z) -le 1e-11 -and [Math]::Abs((Num $p 'area_m2')/$area-1) -le 1e-6 -and @(Need $p 'avg_Ex_re_im').Count -eq 2 -and @(Need $p 'avg_Ey_re_im').Count -eq 2 -and (True (Need $p 'dephased')) -and (Need $p 'phase_convention') -eq 'exp(+i*omega*t); multiply exp(+i*k_parallel dot r_parallel)'}}
  Assert-Vector (Need $m 'incident_k_rad_per_m') $incidence.k_in .1 'physical incident wavevector';$test.incident_direction=$true
  $source=Need $m 'explicit_global_source_binding';$test.explicit_source_binding=(True (Need $source 'verified')) -and (Need $source 'proof_schema') -eq 'ANGLE_BAND_SOURCE_POWER_V2'
  foreach($x in @(@('lambda_um','lambda_um'),@('theta_deg','theta_deg'),@('phi_global_deg','phi_deg'))){Assert-EqualNumber $source $x[0] (Num $case $x[1])}
  Assert-Vector (Need $source 'k_in_requested_rad_per_m') $incidence.k_in .1 'source requested incoming k';Assert-Vector (Need $source 'k_out_requested_rad_per_m') $incidence.k_out .1 'source requested outgoing k';Assert-Vector (Need $source 'incoming_E_unit_global') $incidence.E_in 1e-12 'incoming case-bound p/s polarization';Assert-Vector (Need $source 'outgoing_E_unit_global') $incidence.E_out 1e-12 'outgoing case-bound p/s polarization'
  $requested=@(Need $source 'requested');foreach($key in @('E0_readback','Eampl_readback')){$vals=@(Need $source $key);$test.explicit_source_binding=$test.explicit_source_binding -and $vals.Count -eq 3 -and ($vals -join '|') -eq ($requested -join '|')}
  $test.explicit_source_binding=$test.explicit_source_binding -and $requested.Count -eq 3
  if($requested.Count -eq 3){for($j=0;$j -lt 3;$j++){$v=[string]$requested[$j];if($v -notmatch '^([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)(?:\[V/m\])?$'){$test.explicit_source_binding=$false;continue};$n=[double]::Parse($Matches[1],[Globalization.CultureInfo]::InvariantCulture);$test.explicit_source_binding=$test.explicit_source_binding -and [Math]::Abs($n-$incidence.E_out[$j]) -le 1e-12}}
  $modeAxis=@(Tangential-ModeAxis $case);$modeWant=@(($modeAxis[0]*$modeAxis[0]),($modeAxis[1]*$modeAxis[1]),0.);

  $mode=@(Need $m 'normalized_source_mode_E_component_fractions_xyz');$test.mode_source=((Need $m 'normalized_source_mode_component_scope') -eq 'tangential_port_mode' -and $mode.Count -eq 3 -and (Num $m 'normalized_source_mode_target_component_fraction') -ge .999999 -and (Num $m 'normalized_source_mode_orthogonal_leakage_fraction') -ge -1e-10 -and (Num $m 'normalized_source_mode_orthogonal_leakage_fraction') -le 1e-6)
  Assert-Vector (Need $m 'normalized_source_mode_target_tangential_unit_axis') $modeAxis 1e-12 'case-bound true p/s global-phi tangential mode target axis'
  if(!(True (Need $m 'source_mode_projection_is_global_phi_aware'))){throw 'Native source mode projection does not establish global phi'}
  $componentTolerance=if((Num $case 'phi_deg') -eq 0){.000001}else{.0021};if($mode.Count -eq 3){for($j=0;$j -lt 3;$j++){$test.mode_source=$test.mode_source -and (Number $mode[$j]) -and [Math]::Abs([double]$mode[$j]-$modeWant[$j]) -le $componentTolerance}}
  # Match measured global wavevectors to the independent lattice SET. Port
  # order labels need not use the same reference orientation as this basis.
  $modes=Need $case 'expected_modes';$spatialCount=[int](Num $modes 'spatial_propagating_orders');$channelCount=[int](Num $modes 'propagating_polarization_channels');$allCount=[int](Num $modes 'all_polarization_channels')
  $orders=@(Need $m 'order_powers');$sum=0.;$propCount=0;$keySeen=@{};$ports=@{};$spatial=@{};$orderOK=$orders.Count -eq $allCount;$expected=@(Expected-Lattice $case);$zeroIndex=-1
  for($j=0;$j -lt $expected.Count;$j++){if($expected[$j].m -eq 0 -and $expected[$j].n -eq 0){$zeroIndex=$j}}
  $orderOK=$orderOK -and $expected.Count -eq $spatialCount -and $zeroIndex -ge 0
  foreach($x in $orders){
    $key=[string](Need $x 'key');$pn=[string](Need $x 'port_name');if($keySeen.ContainsKey($key) -or $ports.ContainsKey($pn)){$orderOK=$false};$keySeen[$key]=$true;$ports[$pn]=$true
    if($key -notmatch '^(-?\d+),(-?\d+),(.+)$'){$orderOK=$false;continue};$component=$Matches[3];$index=-1
    if($key -eq '0,0,co'){$orderOK=$orderOK -and $pn -eq '1';$prop=$true;$index=$zeroIndex}else{
      $kx=Num $x 'kx_rad_per_m';$ky=Num $x 'ky_rad_per_m';$kz=Num $x 'real_kz_rad_per_m';$ikz=Num $x 'imag_kz_rad_per_m';$prop=($kx*$kx+$ky*$ky) -lt $k0*$k0*(1-1e-12)
      $actualProp=Need $x 'propagating';$orderOK=$orderOK -and ($actualProp -is [bool]) -and $actualProp -eq $prop
      if($prop){for($j=0;$j -lt $expected.Count;$j++){$dx=$kx-$expected[$j].kx;$dy=$ky-$expected[$j].ky;if([Math]::Sqrt($dx*$dx+$dy*$dy) -lt .1){if($index -ge 0){$orderOK=$false};$index=$j}};$orderOK=$orderOK -and $index -ge 0 -and $kz -gt 0 -and [Math]::Abs($ikz) -lt .1 -and [Math]::Abs(($kz*$kz+$kx*$kx+$ky*$ky)/($k0*$k0)-1) -lt 1e-7}
    }
    if($prop){$power=Num $x 'reflected_power_over_Pin';$orderOK=$orderOK -and $power -ge 0;$sum+=$power;$propCount++;$pair=[string]$index;if(!$spatial.ContainsKey($pair)){$spatial[$pair]=@()};$spatial[$pair]+=$component}
  }
  $test.channels=($orderOK -and $propCount -eq $channelCount -and $spatial.Count -eq $spatialCount);foreach($key in $spatial.Keys){$pair=@($spatial[$key]|Sort-Object);$wanted=if($key -eq [string]$zeroIndex){'co|orthogonal'}else{'inplane|outofplane'};$test.channels=$test.channels -and $pair.Count -eq 2 -and ($pair -join '|') -eq $wanted}
  $test.channel_sum=([Math]::Abs($sum-$r) -le 1e-10)
  $counts=@(Need $m 'actual_independent_spatial_channel_counts');$test.native_actual_channel_counts=($counts.Count -eq $spatialCount -and @($counts|Where-Object {$_ -ne 2}).Count -eq 0)
  $sides=@(Need $m 'side_pair_powers');$side=0.;$test.side_pairs=($sides.Count -eq 3);foreach($x in $sides){$sumSide=(Num $x 'source_W')+(Num $x 'destination_W');$side=[Math]::Max($side,[Math]::Abs($sumSide));$test.side_pairs=$test.side_pairs -and [Math]::Abs($sumSide-(Num $x 'sum_W')) -le 1e-10};$test.side_power=($side -le .002 -and [Math]::Abs($side-(Num $m 'side_pair_power_mismatch_over_Pin')) -le 1e-10)
  $failed=@($test.Keys|Where-Object {!$test[$_]});$nativePassed=$m.state -eq 'PASSED_FLOOR_FILL_CAVITY_DIAGNOSTIC_ONLY';$allNative=@($checks.PSObject.Properties|Where-Object {!(True $_.Value)}).Count -eq 0
  if($m.state -notin @('PASSED_FLOOR_FILL_CAVITY_DIAGNOSTIC_ONLY','FAILED_FLOOR_FILL_CAVITY_DIAGNOSTIC') -or $nativePassed -ne $allNative){throw 'Native result state and boolean check flags disagree'}
  $strict=[ordered]@{mode_vs_net_0p002=($modeGap -le .002);port_vs_loss_0p002=($gap -le .002);net_vs_loss_0p002=($netgap -le .002);Ti_E_vs_H_0p002=($tiGap -le .002);film_manual_vs_Qh_0p002=($filmGap -le .002)}
  $disposition=Exploration-Disposition $test
  return [ordered]@{schema=1;case_id=$CaseId;accepted=($failed.Count -eq 0);state=$(if($failed.Count -eq 0){'EXPLORATORY_SINGLE_POINT_ACCEPTED_NOT_MESH_CONVERGED'}else{'EXPLORATORY_REJECTED'});checks=$test;failed_checks=$failed;strict_original_0p002_diagnostics=$strict;mode_net_abs=$modeGap;port_loss_abs=$gap;port_loss_relative=$relative;net_loss_abs=$netgap;net_loss_relative=$netRelative;Ti_E_vs_H_abs=$tiGap;Ti_E_vs_H_relative=$tiRelative;film_manual_vs_Qh_abs=$filmGap;film_manual_vs_Qh_relative=$filmRelative;raw_metrics_changed=$false;formal_qualification=$false;mesh_convergence_established=$false;bare_C0_is_not_an_accepted_performance_reference=$true;scope='G99 3um bottom-fill single-point optical proxy';exploration_continue_allowed=$disposition.exploration_continue_allowed;disposition=$disposition.disposition;warning_checks=$disposition.warning_checks;hard_failed_checks=$disposition.hard_failed_checks;warning_policy=$disposition.warning_policy}
}
if($ResourceSelfTest -or $UnitSelfTest){
  $tests=@(@(100,100,0,0,0,0,$null),@(129,100,0,0,1,0,$null),@(129,100,1,0,2,0,'OWNED_RSS_HARD_LIMIT'),@(100,63,0,0,0,1,$null),@(100,63,0,1,0,2,'HOST_FREE_MEMORY_RESERVE'),@(128,64,0,0,0,0,$null))
  foreach($x in $tests){$v=Guard ([long]$x[0]*1GB) $x[1] $x[2] $x[3];if($v.hard -ne $x[4] -or $v.low -ne $x[5] -or $v.reason -ne $x[6]){throw 'Resource self-test failed'}}
  if($UnitSelfTest){




    $bindings=@(@('G99_A3_P80_t3_l10p5_th0_phi0_s',10.5,0,0,'s',163,398,98.09681325894577,'m4_air29_wft3',3,34.65468216965557),@('G99_A0facets_P80_t3_l10p5_th0_phi0_s',10.5,0,0,'s',163,398,98.09681325894577,'m4_air29_wft3',3,34.65468216965557))
    foreach($b in $bindings){
     $row=[pscustomobject](Case-Contract $b[0]);if($row.bottom_flat_to_flat_um -ne $b[10] -or $row.height_um -ne $b[7] -or $row.film_um -ne $b[9] -or $row.lambda_um -ne $b[1] -or $row.theta_deg -ne $b[2] -or $row.phi_deg -ne $b[3] -or $row.polarization -ne $b[4] -or $row.mesh -ne $b[8] -or @(Expected-Lattice $row).Count -ne $b[5] -or 2*@(Expected-Lattice $row 1.1).Count -ne $b[6]){throw 'Closed G99 source/lattice contract failed'}
     $inc=Incidence $row;$axis=@(Tangential-ModeAxis $row);$dot=0.;$axisDot=0.;for($j=0;$j -lt 3;$j++){$dot+=$inc.k_in[$j]*$inc.E_in[$j];$axisDot+=$inc.k_in[$j]*$axis[$j]};if([Math]::Abs($dot) -gt 1e-8){throw 'Requested source not transverse'}
     $pr=$row.phi_deg*[Math]::PI/180.;$wantAxis=if($row.polarization -eq 's'){@((-[Math]::Sin($pr)),[Math]::Cos($pr),0.)}else{@([Math]::Cos($pr),[Math]::Sin($pr),0.)};Assert-Vector $axis $wantAxis 1e-12 'independent selected p/s tangential mode axis'
     if($row.polarization -eq 's' -and [Math]::Abs($axisDot) -gt 1e-8){throw 's tangential mode must be orthogonal to in-plane k'}
    }
    $srow=[pscustomobject]@{phi_deg=30;polarization='s'};Assert-Vector (Tangential-ModeAxis $srow) @(-.5,([Math]::Sqrt(3)/2),0.) 1e-12 's phi30 is not the p target axis'
    foreach($id in @('V1_A3H80_t3_l10p5_th0_p_localmesh','V2_A3H80_t3_l13_th15_p','V5_A3H80_t3_l10p5_th30_phi30_p')){$rejected=$false;try{Assert-ExecutableCase $id}catch{$rejected=$true};if(!$rejected){throw 'Prior/unapproved case executable'}}
    [pscustomobject]@{state='PASS_CONTROLLER_UNIT_SELF_TEST';resource_guard_cases=$tests.Count;closed_physical_cases=2;p_cases=0;s_cases=2;independent_shifted_lattice_cases=2;genuine_s_phi30_axis_checked=$true;localmesh_cases=0;old_case_rejections=3;cores_per_case=32;concurrency=1;solver_executed=$false;remote_writes=$false}|ConvertTo-Json;return

  }
  [pscustomobject]@{state='PASS_RESOURCE_GUARD_SELF_TEST';cases=$tests.Count;solver_executed=$false;remote_writes=$false}|ConvertTo-Json;return
}
$stage=[IO.Path]::GetFullPath($StageRoot)
if($stage -ne '__STAGE_ROOT__'){throw 'Unexpected isolated floor-fill angular/band stage root'}
if($PlanSelfTest){Assert-Manifest;$plan=Json (Join-Path $stage 'plan.json');Assert-Plan;[pscustomobject]@{state='PASS_FROZEN_CASE_PLAN_SELF_TEST';case_id=$CaseId;polarization=$case.polarization;solver_executed=$false;remote_writes=$false;window_not_required_for_plan_selftest=$true}|ConvertTo-Json;return}
if($EvaluationSelfTest){
  Assert-Manifest;$plan=Json (Join-Path $stage 'plan.json');Assert-Plan
  $fixture=Json (Join-Path $stage 'warn_fixture_metrics.json');$checks=[ordered]@{}
  foreach($p in (Need $fixture 'checks').PSObject.Properties){$checks['native_'+$p.Name]=True $p.Value}
  $net=Num $fixture 'A_flux';$loss=Num $fixture 'A_loss_total';$r=Num $fixture 'R_total';$top=Num $fixture 'Pout_top_net_W'
  $checks.mode_vs_net=([Math]::Abs($r-1-$top) -le .005);$checks.net_vs_loss_abs=([Math]::Abs($net-$loss) -le .005);$checks.net_vs_loss_relative=([Math]::Abs($net-$loss)/[Math]::Max([Math]::Abs($net),1e-12) -le .02)
  $checks.material_permittivity=$true;$checks.incident_power=$true
  $actual=Exploration-Disposition $checks
  if($actual.accepted -or !$actual.exploration_continue_allowed -or $actual.disposition -ne 'COMPLETE_WITH_FLUX_WARNING' -or @($actual.warning_checks).Count -ne 6 -or @($actual.hard_failed_checks).Count -ne 0){throw 'Real-Q0 fixture must retain accepted=false and exactly six known nonblocking flux failures'}
  $results=[Collections.Generic.List[object]]::new();$results.Add([ordered]@{fixture='REAL_SAVED_Q0_METRICS';classification=$actual})
  foreach($hardName in @('incident_power','material_permittivity','native_Ti_H_vs_E_loss_relative','native_periodic_side_power','unknown_future_check')){
    $injected=[ordered]@{};foreach($key in $checks.Keys){$injected[$key]=$checks[$key]};$injected[$hardName]=$false;$d=Exploration-Disposition $injected
    if($d.accepted -or $d.exploration_continue_allowed -or $d.disposition -ne 'STOP_HARD_FAILURE' -or $hardName -notin $d.hard_failed_checks){throw ('Hard failure was incorrectly treated as warning: '+$hardName)}
    $results.Add([ordered]@{injected_hard_failure=$hardName;classification=$d})
  }
  $clean=[ordered]@{};foreach($key in $checks.Keys){$clean[$key]=$true};$d=Exploration-Disposition $clean;if(!$d.accepted -or !$d.exploration_continue_allowed -or @($d.failed_checks).Count){throw 'Clean fixture acceptance definition changed'}
  [ordered]@{state='PASS_EVALUATION_WARNING_POLICY_SELF_TEST';real_fixture_sha256=Hash (Join-Path $stage 'warn_fixture_metrics.json');results=$results.ToArray();all_checks_true_accepted=$true;accepted_definition_unchanged=$true;solver_executed=$false;remote_writes=$false}|ConvertTo-Json -Depth 20;return
}
if($ExpectedManifestSha -notmatch '^[a-fA-F0-9]{64}$' -or $CaseId -notmatch '^[A-Za-z0-9][A-Za-z0-9_-]{0,79}$' -or $AttemptId -notmatch '^[A-Za-z0-9][A-Za-z0-9_-]{0,47}$'){throw 'Unsafe identifiers or missing expected manifest SHA256'}
Assert-ExecutableCase $CaseId
$short=if($Mode -eq 'MeshOnly'){'mesh'}else{'solve'};$suffix=if($AttemptId -eq 'initial'){''}else{'_'+$AttemptId}
$statusPath=Join-Path $stage ('controller_'+$CaseId+'_'+$short+$suffix+'.json')
if(Test-Path -LiteralPath $statusPath){throw 'Attempt status exists; no implicit restart'}
$lock=[IO.File]::Open((Join-Path $stage 'controller.lock'),[IO.FileMode]::OpenOrCreate,[IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
try{
  if(Test-Path -LiteralPath $statusPath){throw 'Attempt appeared while acquiring exclusive stage lock'}
  Initialize-ExecutionContext;Save-Status 'STARTING'
  $reviewed=$null;if($Mode -eq 'SolveApproved'){$reviewed=Assert-Admission}
  $r=Wait-Resources;if($Mode -eq 'SolveApproved'){$reviewed=Assert-Admission};$parent=Join-Path $stage $CaseId;if(!(Test-Path -LiteralPath $parent)){New-Item -ItemType Directory -Path $parent -ErrorAction Stop|Out-Null}
  $leaf=if($Mode -eq 'MeshOnly'){'mesh_preflight'+$suffix}else{'solve'+$suffix};$dir=Join-Path $parent $leaf
  if(Test-Path -LiteralPath $dir){throw 'Attempt directory exists; refusing overwrite or retry'};New-Item -ItemType Directory -Path $dir -ErrorAction Stop|Out-Null
  $receipt=Write-CaseEntry $dir ($Mode -eq 'MeshOnly') $reviewed
  Assert-Window
  $active=[pscustomobject]@{case_id=$CaseId;directory=$dir;process_id=$null;owned_rss_bytes=0L;owned_process_ids=@();hard_samples=0;host_low_samples=0;native_progress=$null};Save-Status 'RUNNING' $r
  $exitCode=Run-OwnedCmd $dir (Join-Path $dir 'run_case.cmd');$active=$null
  $compileExit=[int]([IO.File]::ReadAllText((Join-Path $dir 'compile_exitcode.txt'),[Text.Encoding]::ASCII).Trim())
  New-Json (Join-Path $dir 'native_exit_receipt.json') ([ordered]@{case_id=$CaseId;compile_exit_code=$compileExit;native_exit_code=$exitCode;compile_log_tail_16KiB=Log-Tail (Join-Path $dir 'compile.stdout.log');native_log_tail_16KiB=Log-Tail (Join-Path $dir 'stdout.log');logs_retained=$true})
  if($compileExit -ne 0){throw ('Native compilation failed; compile exitcode '+$compileExit)}
  if((Hash (Join-Path $dir 'CaseEntry.java')) -ne $receipt.generated_source_sha256){throw 'Generated source changed'}
  $compiled=[IO.File]::ReadAllText((Join-Path $dir 'compiled_case.sha256'),[Text.Encoding]::ASCII).Trim();if($compiled -notmatch '^[a-f0-9]{64}$' -or (Hash (Join-Path $dir 'CaseEntry.class')) -ne $compiled){throw 'Compiled class fingerprint missing or changed'}
  $receipt.compiled_case_sha256=$compiled;Write-Json (Join-Path $dir 'case_receipt.json') $receipt;Assert-Manifest
  if($Mode -eq 'SolveApproved' -and (Hash (Admission-Path)) -ne $admissionHash){throw 'Admission changed during physical attempt'}
  $stdout=Join-Path $dir 'stdout.log';$nativeFailure=Native-Failure $stdout;if($null -ne $nativeFailure){New-Json (Join-Path $dir 'native_failure.json') $nativeFailure;throw (Native-FailureMessage $nativeFailure)}
  if($exitCode -ne 0){throw ('Native process failed; exitcode '+$exitCode+'; complete native stdout retained')}
  $geometry=Marker $stdout 'CODEX_GEOMETRY_JSON=';$mesh=Marker $stdout 'CODEX_MESH_JSON='
  New-Json (Join-Path $dir 'geometry.json') $geometry;New-Json (Join-Path $dir 'mesh.json') $mesh
  $portProof=Marker $stdout 'CODEX_PRESOLVE_PORT_BASIS_JSON=';Assert-PresolvePortProof $portProof
  if(($portProof|ConvertTo-Json -Depth 30 -Compress) -ne ((Need $mesh 'presolve_port_wavevector_proof')|ConvertTo-Json -Depth 30 -Compress)){throw 'Independent presolve-port marker and saved mesh inventory disagree'}
  New-Json (Join-Path $dir 'presolve_port_wavevectors.json') $portProof;$receipt['presolve_port_wavevector_proof_sha256']=Hash (Join-Path $dir 'presolve_port_wavevectors.json');Write-Json (Join-Path $dir 'case_receipt.json') $receipt
  $ooc=Marker $stdout 'CODEX_LINEAR_SOLVER_OVERRIDE_JSON=';Assert-OOC $ooc
  if(($ooc|ConvertTo-Json -Depth 20 -Compress) -ne ($mesh.linear_solver_override|ConvertTo-Json -Depth 20 -Compress)){throw 'Independent OOC marker and saved mesh inventory disagree'}
  New-Json (Join-Path $dir 'linear_solver_override.json') $ooc;$receipt.linear_solver_override_sha256=Hash (Join-Path $dir 'linear_solver_override.json');Write-Json (Join-Path $dir 'case_receipt.json') $receipt;Assert-Mesh $mesh $geometry
  if(!(Test-Path -LiteralPath (Join-Path $dir 'mesh_only.mph') -PathType Leaf)){throw 'No persisted mesh_only.mph'}
  if($Mode -eq 'MeshOnly'){
    $done=Marker $stdout 'CODEX_DONE_JSON=';if($done.case_id -ne $CaseId -or !(True (Need $done 'mesh_only')) -or (Need $done 'solver_executed') -isnot [bool] -or $done.solver_executed -or (Need $done 'model_file') -ne (Join-Path $dir 'mesh_only.mph')){throw 'Mesh DONE evidence missing or inconsistent'}
    New-Json (Join-Path $dir 'done.json') $done;$completed=[ordered]@{case_id=$CaseId;directory=$dir;exit_code=$exitCode;eligible_for_manual_review=$true;degrees_of_freedom=$mesh.degrees_of_freedom;tetrahedra=$mesh.tetrahedra;mesh_model_sha256=Hash (Join-Path $dir 'mesh_only.mph');solver_executed=$false};Save-Status 'COMPLETE_MESH_PREFLIGHT_AWAITING_EXPLICIT_ADMISSION'
  }else{
    if($mesh.degrees_of_freedom -ne $reviewed.degrees_of_freedom -or $mesh.tetrahedra -ne $reviewed.tetrahedra -or $mesh.mesh_fingerprint_sha256 -ne $reviewed.mesh_fingerprint_sha256){throw 'Solve mesh differs from explicitly reviewed actual inventory or exact fingerprint'}
    $m=Marker $stdout 'CODEX_METRICS_JSON=';New-Json (Join-Path $dir 'metrics.json') $m
    if(($ooc|ConvertTo-Json -Depth 20 -Compress) -ne ($m.linear_solver_override|ConvertTo-Json -Depth 20 -Compress)){throw 'Physical metrics OOC readback differs from the actual admitted configuration'}
    $solved=Marker $stdout 'CODEX_SOLVED_ARTIFACT_JSON=';$done=Marker $stdout 'CODEX_DONE_JSON='
    if($solved.case_id -ne $CaseId -or !(True (Need $solved 'solver_executed')) -or (Need $solved 'model_path') -ne (Join-Path $dir 'solved_unevaluated.mph') -or $done.case_id -ne $CaseId -or (Need $done 'mesh_only') -isnot [bool] -or $done.mesh_only -or !(True (Need $done 'solver_executed')) -or (Need $done 'model_file') -ne (Join-Path $dir 'model.mph')){throw 'Solved/checkpoint/DONE identities do not agree'}
    foreach($n in @('solved_unevaluated.mph','model.mph')){if(!(Test-Path -LiteralPath (Join-Path $dir $n) -PathType Leaf)){throw ('Missing saved physical artifact '+$n)}}
    New-Json (Join-Path $dir 'solved_artifact.json') $solved;New-Json (Join-Path $dir 'done.json') $done;$evaluation=Evaluate-Metrics $m $mesh $geometry;New-Json (Join-Path $dir 'evaluation.json') $evaluation
    $artifactHashes=[ordered]@{case_id=$CaseId;model_sha256=Hash (Join-Path $dir 'model.mph');solved_checkpoint_sha256=Hash (Join-Path $dir 'solved_unevaluated.mph');metrics_sha256=Hash (Join-Path $dir 'metrics.json');evaluation_sha256=Hash (Join-Path $dir 'evaluation.json');geometry_sha256=Hash (Join-Path $dir 'geometry.json');mesh_sha256=Hash (Join-Path $dir 'mesh.json');receipt_sha256=Hash (Join-Path $dir 'case_receipt.json')}
    New-Json (Join-Path $dir 'physical_artifact_hashes.json') $artifactHashes
    if(!$evaluation.exploration_continue_allowed){throw 'Floor-fill cavity hard validation failed; raw metrics, original flags, evaluation, saved model fingerprints and all failure evidence retained'}
    $completed=[ordered]@{case_id=$CaseId;directory=$dir;exit_code=$exitCode;metrics_sha256=$artifactHashes.metrics_sha256;model_sha256=$artifactHashes.model_sha256;solved_checkpoint_sha256=$artifactHashes.solved_checkpoint_sha256;evaluation_sha256=$artifactHashes.evaluation_sha256;physical_artifact_hashes_sha256=Hash (Join-Path $dir 'physical_artifact_hashes.json');evaluation=$evaluation;solver_executed=$true;formal_qualification=$false;mesh_convergence_established=$false}
    if($evaluation.accepted){Save-Status 'COMPLETE_FLOOR_FILL_CAVITY_DIAGNOSTIC_ONLY'}else{Save-Status 'COMPLETE_FLOOR_FILL_WITH_FLUX_WARNING'}
  }
}catch{
  $message=$_.Exception.Message;try{Stop-OwnedIdentities}catch{$message+='; OWNED_CLEANUP_REVIEW: '+$_.Exception.Message};$active=$null
  try{Save-Status 'STOPPED_REVIEW_REQUIRED' $null $message}catch{};throw
}finally{if($null -ne $lock){$lock.Dispose()}}
