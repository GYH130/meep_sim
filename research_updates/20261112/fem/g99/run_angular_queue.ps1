param([Parameter(Mandatory=$true)][string]$ExpectedManifestSha,[string]$StageRoot='__STAGE_ROOT__',[switch]$QueueSelfTest,[switch]$PrerequisiteSelfTest)
$ErrorActionPreference='Stop';$OutputEncoding=[Console]::OutputEncoding=[Text.UTF8Encoding]::new($false)
$stage=[IO.Path]::GetFullPath($StageRoot);$FixedStage='__STAGE_ROOT__'
if ($stage -ne $FixedStage -or $ExpectedManifestSha -notmatch '^[a-fA-F0-9]{64}$' -or
    [IO.Path]::GetFullPath($PSCommandPath) -ne (Join-Path $stage 'run_angular_queue.ps1')) {
    throw 'Frozen isolated G99 stage/SHA required'
}
$helperPath=Join-Path $stage 'queue_review_helpers.ps1'
$manifestPath=Join-Path $stage 'launch_manifest.json'
if ((Get-FileHash -LiteralPath $manifestPath -Algorithm SHA256).Hash.ToLowerInvariant() -ne $ExpectedManifestSha.ToLowerInvariant()) {
    throw 'Frozen G99 manifest differs'
}
$initialManifest=[IO.File]::ReadAllText($manifestPath,[Text.Encoding]::UTF8) | ConvertFrom-Json
if ((Get-FileHash -LiteralPath $helperPath -Algorithm SHA256).Hash.ToLowerInvariant() -ne
    $initialManifest.source_hashes.'queue_review_helpers.ps1') {
    throw 'Frozen G99 review helpers differ'
}
. $helperPath
$CaseIds=@('G99_A3_P80_t3_l10p5_th0_phi0_s','G99_A0facets_P80_t3_l10p5_th0_phi0_s');$MainCaseIds=@($CaseIds);$ConditionalCaseIds=@()
$AllowedWarnings=@('native_R_orders_vs_net_flux','native_flux_vs_total_loss_absolute','native_flux_vs_total_loss_relative','mode_vs_net','net_vs_loss_abs','net_vs_loss_relative')
$manifestHash=$null;$manifest=$null;$windowHash=$null;$window=$null;$lock=$null;$statusOwned=$false;$startedUtc=$null;$finishedUtc=$null;$deadlineUtc=$null
$statusPath=Join-Path $stage 'queue_status.json';$windowPath=Join-Path $stage 'window.json';$controller=Join-Path $stage 'run_floorfilm_controller.ps1'
$completed=[Collections.Generic.List[object]]::new();$skipped=[Collections.Generic.List[object]]::new();$skippedCases=@{};$active=$null;$peakRss=0L;$pairComparison=$null;$plan=$null;$meshByThickness=@{};$meshByHeight=@{};$preflight=[Collections.Generic.List[object]]::new();$admissions=@{};$pairAdmissionHash=$null;$resourceSnapshot=$null;$calibrationReview=$null;$diagnosticReview=$null;$shortwaveReview=$null;$CurrentCaseId=$CaseIds[0]
# Stage-local override of the pinned parent helper: same body, authorized24GiB cache.
function Assert-Ooc($Evidence,[string]$CaseId){
  if((Need $Evidence 'case_id') -ne $CaseId -or (Need $Evidence 'state') -ne 'VERIFIED_OOC_CONFIGURATION_BEFORE_PHYSICS' -or (Need $Evidence 'method') -ne 'EXACT_DIRECT_OUT_OF_CORE_NO_BLR' -or (Need $Evidence 'stationary_tag') -ne 's1'){throw 'Missing case-bound exact OOC evidence'}
  foreach($n in @('readback_verified','active_selection_preserved','oocmemory_is_not_process_RSS_cap')){if(!(Is-True (Need $Evidence $n))){throw ('OOC verification not true: '+$n)}}
  if(!(Is-False (Need $Evidence 'physical_thresholds_changed'))){throw 'OOC evidence reports changed physical thresholds'}
  foreach($n in @('active_selection_before','active_selection_after')){$a=Need $Evidence $n;if(@($a.PSObject.Properties).Count -ne 1 -or (Need $a 'fc1') -ne 'd1'){throw 'Original active fc1->d1 selection changed'}}
  Equal-Number $Evidence 'direct_feature_count' 2
  $changes=@(Need $Evidence 'direct_changes');if($changes.Count -ne 2){throw 'Exactly two original Direct readbacks required'}
  foreach($pair in @(@('dDef','mumps'),@('d1','pardiso'))){
    $rows=@($changes|Where-Object {$_.tag -eq $pair[0]})
    if($rows.Count -ne 1 -or !(Is-True (Need $rows[0] 'readback_verified'))){throw 'Missing/duplicate Direct feature readback'}
    $before=Need $rows[0] 'before';$after=Need $rows[0] 'after'
    if((Need $before 'linsolver') -ne $pair[1] -or (Need $after 'linsolver') -ne $pair[1] -or (Need $after 'ooc') -ne 'on' -or (Need $after 'incore') -ne 'manual' -or (Need $after 'mumpsblr') -ne 'off'){throw 'Direct OOC algorithm/settings changed'}
    Equal-Number $after 'oocmemory_MB' 24576
  }
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
function Mesh-FamilyKey($Row){return [string](Need $Row 'source_geometry_id')}

function Controller-Incidence($Row){
 $theta=(Num $Row 'theta_deg')*[Math]::PI/180.;$phi=(Num $Row 'phi_deg')*[Math]::PI/180.;$k0=2*[Math]::PI/((Num $Row 'lambda_um')*1e-6)
 $st=[Math]::Sin($theta);$ct=[Math]::Cos($theta);$sp=[Math]::Sin($phi);$cp=[Math]::Cos($phi);$pol=[string](Need $Row 'polarization')
 if($pol -eq 'p'){$ein=@(($ct*$cp),($ct*$sp),$st);$eout=@(($ct*$cp),($ct*$sp),(-$st))}elseif($pol -eq 's'){$ein=@((-$sp),$cp,0.);$eout=@((-$sp),$cp,0.)}else{throw 'Only explicit p/s polarization allowed'}
 return [ordered]@{k0=$k0;sin_theta=$st;cos_theta=$ct;k_in=@(($k0*$st*$cp),($k0*$st*$sp),(-$k0*$ct));k_out=@(($k0*$st*$cp),($k0*$st*$sp),($k0*$ct));E_in=$ein;E_out=$eout}
}
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

function Assert-Manifest {
 $p=Join-Path $stage 'launch_manifest.json';$h=Hash $p;if($h -ne $ExpectedManifestSha.ToLowerInvariant() -or ($null -ne $script:manifestHash -and $h -ne $script:manifestHash)){throw 'Frozen manifest changed'};$m=Json $p
 if((Need $m 'purpose') -ne 'G99_P80_MATCHED_WALL_FULLWAVE_PILOT' -or (Need $m 'stage_root') -ne $stage -or (@(Need $m 'physical_case_ids') -join '|') -ne ($CaseIds -join '|') -or @(Need $m 'mesh_only_case_ids').Count -ne 0 -or (@(Need $m 'main_case_ids') -join '|') -ne ($MainCaseIds -join '|') -or (@(Need $m 'conditional_case_ids') -join '|') -ne ($ConditionalCaseIds -join '|')){throw 'Manifest scope differs'}
 if(!(Is-True (Need $m 'physical_thresholds_unchanged'))){throw 'Scientific thresholds changed'};Equal-Number $m 'automatic_retries' 0
 $s=Need $m 'source_hashes';foreach($n in @('CavityFloorFilm.java','run_floorfilm_controller.ps1','run_angular_queue.ps1','run_flat_controls.ps1','FlatObliqueAdmission.java','plan.json','material.json','warn_fixture_metrics.json','source_lineage_audit.json','queue_review_helpers.ps1')){[void](Need $s $n)}
 foreach($p in $s.PSObject.Properties){if([IO.Path]::GetFileName($p.Name) -ne $p.Name -or $p.Name -match '[\\/]' -or $p.Value -notmatch '^[a-f0-9]{64}$' -or (Hash (Join-Path $stage $p.Name)) -ne $p.Value){throw ('Frozen source changed: '+$p.Name)}}
 $script:manifestHash=$h;$script:manifest=$m
}

function Independent-OrderSet([double]$Wavelength,[double]$Theta,[double]$FrequencyRatio,[double]$Phi=0.){
 $k=2*[Math]::PI/$Wavelength*$FrequencyRatio;$tr=$Theta*[Math]::PI/180.;$pr=$Phi*[Math]::PI/180.;$kx0=$k*[Math]::Sin($tr)*[Math]::Cos($pr);$ky0=$k*[Math]::Sin($tr)*[Math]::Sin($pr);$bx=2*[Math]::PI/80.;$by=$bx/[Math]::Sqrt(3);$keys=[Collections.Generic.List[string]]::new();$margin=[double]::PositiveInfinity
 for($m=-40;$m -le 40;$m++){for($n=-40;$n -le 40;$n++){$x=$kx0+$m*$bx;$y=$ky0+(-$m+2*$n)*$by;$v=($x*$x+$y*$y)/($k*$k);$margin=[Math]::Min($margin,[Math]::Abs(1-$v));if($v -lt 1-1e-12){if([Math]::Abs($m) -eq 40 -or [Math]::Abs($n) -eq 40){throw 'Independent lattice reached enumeration boundary'};$keys.Add(([string]$m+','+[string]$n))}}};if($margin -lt 1e-8){throw 'Near-grazing channel requires separate review'};return [pscustomobject]@{keys=@($keys.ToArray());minimum_grazing_q_distance=$margin}
}
function Assert-IndependentModes($Row){
  $e=Need $Row 'expected_modes';$physical=Independent-OrderSet (Num $Row 'lambda_um') (Num $Row 'theta_deg') 1. (Num $Row 'phi_deg');$padding=Independent-OrderSet (Num $Row 'lambda_um') (Num $Row 'theta_deg') 1.1 (Num $Row 'phi_deg')
  Equal-Number $e 'spatial_propagating_orders' @($physical.keys).Count;Equal-Number $e 'propagating_polarization_channels' (2*@($physical.keys).Count);Equal-Number $e 'all_polarization_channels' (2*@($padding.keys).Count);Equal-Number $e 'padding_generation_frequency_ratio' 1.1
  foreach($pair in @(@('physical_order_indices',$physical),@('padding_order_indices',$padding))){$seen=[Collections.Generic.List[string]]::new();foreach($index in @(Need $e $pair[0])){if(@($index).Count -ne 2){throw 'Invalid frozen independent order index'};$seen.Add(([string]$index[0]+','+[string]$index[1]))};if(($seen.ToArray() -join '|') -ne ($pair[1].keys -join '|')){throw 'Frozen independent reciprocal order list differs'}}
}
function Assert-PresolveProof($Proof,$Row){
 if((Need $Proof 'case_id') -ne (Need $Row 'id') -or (Need $Proof 'proof_schema') -ne 'ANGLE_BAND_SOURCE_POWER_V2' -or (Need $Proof 'state') -ne 'PASS_PRESOLVE_COMPLETE_PORT_BASIS' -or !(Is-False (Need $Proof 'solver_executed')) -or !(Is-False (Need $Proof 'native_field_solver_run_before_proof')) -or (Need $Proof 'source_frequency_expression') -ne 'f0'){throw 'Complete variables-only port proof required before physical solver'}
 foreach($n in @('physical_frequency_restored','incident_wavevector_matches_requested','actual_wavevectors_match_independent_lattice','physical_index_set_exact','padding_index_set_exact','port_names_unique')){if(!(Is-True (Need $Proof $n))){throw ('Presolve port proof flag false: '+$n)}}
 foreach($p in @(@('wavelength_um','lambda_um'),@('theta_deg','theta_deg'),@('phi_global_deg','phi_deg'))){Equal-Number $Proof $p[0] (Num $Row $p[1])}
 Equal-Number $Proof 'generator_superset_frequency_ratio' 1.2;Equal-Number $Proof 'retained_padding_frequency_ratio' 1.1;Equal-Number $Proof 'vector_absolute_tolerance_rad_per_m' .1;Equal-Number $Proof 'dispersion_relative_tolerance' 1e-7
 $physical=Independent-OrderSet (Num $Row 'lambda_um') (Num $Row 'theta_deg') 1. (Num $Row 'phi_deg');$padding=Independent-OrderSet (Num $Row 'lambda_um') (Num $Row 'theta_deg') 1.1 (Num $Row 'phi_deg')
 Equal-Number $Proof 'physical_spatial_orders' @($physical.keys).Count;Equal-Number $Proof 'physical_polarization_channels' (2*@($physical.keys).Count);Equal-Number $Proof 'total_port_channels' (2*@($padding.keys).Count);Equal-Number $Proof 'evanescent_channels' (2*(@($padding.keys).Count-@($physical.keys).Count))
 $k0=2*[Math]::PI/((Num $Row 'lambda_um')*1e-6);$theta=(Num $Row 'theta_deg')*[Math]::PI/180;$want=@((Controller-Incidence $Row).k_in);$k=@(Need $Proof 'actual_incident_k_rad_per_m');$requested=@(Need $Proof 'requested_k_in_rad_per_m')
 if($k.Count -ne 3 -or $requested.Count -ne 3){throw 'Presolve incident vector missing'};$norm=0.;for($j=0;$j -lt 3;$j++){if(!(Number $k[$j]) -or !(Number $requested[$j]) -or [Math]::Abs([double]$k[$j]-$want[$j]) -gt .1 -or [Math]::Abs([double]$requested[$j]-$want[$j]) -gt 1e-7){throw 'Presolve source vector differs from original requested angle'};$norm+=[double]$k[$j]*[double]$k[$j]};if([Math]::Abs([Math]::Sqrt($norm)-$k0) -gt .1 -or $k[2] -ge 0){throw 'Presolve incoming vector magnitude/direction wrong'}
 $rows=@(Need $Proof 'channels');if($rows.Count -ne 2*@($padding.keys).Count){throw 'Presolve actual channel count differs from exact unchanged padding'};$seen=@{};$names=@{};$counts=@{};$pols=@{};$physicalCounts=@{};$bx=2*[Math]::PI/(80e-6);$by=$bx/[Math]::Sqrt(3)
 foreach($r in $rows){
  $cm=Num $r 'canonical_m';$cn=Num $r 'canonical_n';if($cm -ne [Math]::Floor($cm) -or $cn -ne [Math]::Floor($cn) -or (Num $r 'native_m') -ne -$cm -or (Num $r 'native_n') -ne -$cn){throw 'Port reciprocal label mapping wrong'};$key=([string][int]$cm+','+[string][int]$cn)
  $name=[string](Need $r 'port_name');$component=[string](Need $r 'components');if($key -notin $padding.keys -or $names.ContainsKey($name) -or $seen.ContainsKey($key+','+$component)){throw 'Unexpected/duplicate port channel'};$names[$name]=$true;$seen[$key+','+$component]=$true
  $ex=$want[0]+$cm*$bx;$ey=$want[1]+(-$cm+2*$cn)*$by;$kx=Num $r 'kx_rad_per_m';$ky=Num $r 'ky_rad_per_m';$kz=Num $r 'real_kz_rad_per_m';$ikz=Num $r 'imag_kz_rad_per_m';$q=($kx*$kx+$ky*$ky)/($k0*$k0);$qExpected=($ex*$ex+$ey*$ey)/($k0*$k0);$prop=$key -in $physical.keys
  if([Math]::Abs($kx-$ex) -gt .1 -or [Math]::Abs($ky-$ey) -gt .1 -or [Math]::Abs((Num $r 'expected_kx_rad_per_m')-$ex) -gt 1e-7 -or [Math]::Abs((Num $r 'expected_ky_rad_per_m')-$ey) -gt 1e-7 -or [Math]::Abs(($kx*$kx+$ky*$ky+$kz*$kz-$ikz*$ikz)/($k0*$k0)-1) -gt 1e-7 -or !(Is-True (Need $r 'wavevector_verified')) -or (Need $r 'propagating') -isnot [bool] -or [bool]$r.propagating -ne $prop -or (Need $r 'actual_wavevector_propagating') -isnot [bool] -or [bool]$r.actual_wavevector_propagating -ne $prop){throw 'Actual presolve channel wavevector/classification failed'}
  if($prop){if($q -ge 1-1e-10 -or [Math]::Abs($ikz) -gt .1 -or (Num $r 'independent_expected_spatial_index') -ne [Array]::IndexOf([string[]]$physical.keys,$key)){throw 'Propagating mode changed/evanescent or independent index wrong'};$physicalCounts[$key]=1+[int]$physicalCounts[$key]}
  elseif($q -le 1+1e-10 -or $qExpected -le 1+1e-10 -or [Math]::Abs($kz) -gt .1 -or [Math]::Abs([Math]::Abs($ikz)-$k0*[Math]::Sqrt($qExpected-1)) -gt .1 -or (Num $r 'independent_expected_spatial_index') -ne -1){throw 'Padding evanescent mode does not satisfy exact independent vacuum dispersion'}
  $counts[$key]=1+[int]$counts[$key];if(!$pols.ContainsKey($key)){$pols[$key]=[Collections.Generic.List[string]]::new()};$pols[$key].Add($component)
 }
 foreach($key in $padding.keys){$expectedPol=if($key -eq '0,0'){'co|orthogonal'}else{'inplane|outofplane'};if($counts[$key] -ne 2 -or (@($pols[$key]|Sort-Object) -join '|') -ne $expectedPol){throw 'Every original padding direction requires both polarizations'}}
 foreach($key in $physical.keys){if($physicalCounts[$key] -ne 2){throw 'Physical direction missing either polarization'}}
 $countsReported=@(Need $Proof 'actual_independent_spatial_channel_counts');if($countsReported.Count -ne @($physical.keys).Count -or @($countsReported|Where-Object {$_ -ne 2}).Count){throw 'Presolve physical pair count inventory failed'}
 return [ordered]@{state='PASS_PARENT_INDEPENDENT_PRESOLVE_PORT_REVIEW';case_id=$Row.id;physical_spatial_orders=@($physical.keys).Count;physical_polarization_channels=2*@($physical.keys).Count;total_port_channels=$rows.Count;actual_wavevectors_recomputed=$true;unchanged_original_padding=$true;physical_solver_executed=$false}
}
function Assert-WavevectorFitConsistency($Fit,$K,[double]$Lambda,[double]$Theta,[double]$Phi=0.){
  $k0=2*[Math]::PI/($Lambda*1e-6);$thetaRad=$Theta*[Math]::PI/180;$phiRad=$Phi*[Math]::PI/180.;$want=@(($k0*[Math]::Sin($thetaRad)*[Math]::Cos($phiRad)),($k0*[Math]::Sin($thetaRad)*[Math]::Sin($phiRad)),(-$k0*[Math]::Cos($thetaRad)))
  if(@($K).Count -ne 3){throw 'Incident wavevector missing'}
  $norm2=0.;for($j=0;$j -lt 3;$j++){if(!(Number $K[$j]) -or [Math]::Abs([double]$K[$j]-$want[$j]) -gt .1){throw 'Actual native incident global wavevector differs from intended angle'};$norm2+=[double]$K[$j]*[double]$K[$j]}
  $norm=[Math]::Sqrt($norm2);if($norm -le 0 -or [Math]::Abs($norm-$k0) -gt .1){throw 'Actual native vacuum wavevector magnitude differs from requested frequency'}
  $cos=-[double]$K[2]/$norm
  if($cos -le 0 -or $cos -gt 1){throw 'Incoming direction cosine invalid'}
  if(!(Is-True (Need $Fit 'direction_normalized_to_exact_vacuum_k0'))){throw 'Fit must project actual direction to the exact vacuum wavevector magnitude'}
  $raw=@(Need $Fit 'actual_k_in_rad_per_m');$projected=@(Need $Fit 'projection_k_in_rad_per_m')
  if($raw.Count -ne 3 -or $projected.Count -ne 3){throw 'Source-fit actual/projected wavevectors missing'}
  for($j=0;$j -lt 3;$j++){if(!(Number $raw[$j]) -or !(Number $projected[$j]) -or [Math]::Abs([double]$raw[$j]-[double]$K[$j]) -gt 1e-8 -or [Math]::Abs([double]$projected[$j]-[double]$K[$j]*$k0/$norm) -gt 1e-6){throw 'Source-fit raw/projected direction metadata inconsistent'}}
  if([Math]::Abs((Num $Fit 'actual_k_magnitude_rad_per_m')/$norm-1) -gt 1e-12){throw 'Source-fit actual magnitude metadata inconsistent'}
  # Same normalization as the unchanged production Java fit. The component and
  # magnitude checks above independently retain the requested source contract.
  Equal-Number $Fit 'cos_theta' $cos;Equal-Number $Fit 'power_admittance_p_over_vacuum' (1/$cos);Equal-Number $Fit 'power_admittance_s_over_vacuum' $cos
  if([Math]::Abs((Num $Fit 'k0_fixed_rad_per_m')/$k0-1) -gt 1e-12 -or [Math]::Abs((Num $Fit 'kz_fixed_rad_per_m')/($k0*$cos)-1) -gt 1e-12){throw 'Fixed source-fit k0/kz inconsistent with actual direction'}
  return [ordered]@{actual_norm_rad_per_m=$norm;requested_k0_rad_per_m=$k0;direction_cosine=$cos;component_tolerance_rad_per_m=.1;magnitude_tolerance_rad_per_m=.1;normalization='minus_kz_over_vector_norm'}
}
function Assert-ActualSource($Metrics,$Row){
  $proof='ANGLE_BAND_SOURCE_POWER_V2';$fit=Need $Metrics 'internal_E_power_fit';$method='FIXED_KZ_BLOCH_DEPHASED_FIVE_PLANE_TANGENTIAL_SP_COMPLEX_LS_V2'
  if((Need $Metrics 'proof_schema') -ne $proof -or !(Is-True (Need $Metrics 'requested_incidence_verified')) -or !(Is-True (Need $Metrics 'global_phi_verified')) -or !(Is-True (Need $Metrics 'side_floquet_source_verified'))){throw 'Actual native incidence/globalphi/Floquet source proof missing'}
  Equal-Number $Metrics 'port_azimuth_offset_deg' -150
  if((Need $fit 'proof_schema') -ne $proof -or (Need $fit 'method') -ne $method -or (Need $Metrics 'incident_power_method') -ne $method){throw 'Physical source uses stale normal-only power estimator'}
  [void](Assert-WavevectorFitConsistency $fit (Need $Metrics 'incident_k_rad_per_m') (Num $Row 'lambda_um') (Num $Row 'theta_deg') (Num $Row 'phi_deg'))
  foreach($pair in @(@('nfilm','n_film'),@('kfilm','k_film'),@('nti','n_ti'),@('kti','k_ti'))){Equal-Number $Metrics $pair[0] (Num $Row $pair[1])}
  if([Math]::Abs((Num $fit 'Pin_W')-1) -gt .005 -or [Math]::Abs((Num $Metrics 'Pin_independent_internal_E_fit_W')-(Num $fit 'Pin_W')) -gt 1e-10 -or (Num $fit 'maximum_full_or_LOO_design_condition_2') -gt 10 -or (Num $fit 'relative_L2_E_fit_residual') -gt .005 -or (Num $fit 'LOO_max_Pin_relative_error_vs_specified_1W') -gt .005 -or !(Is-True (Need $fit 'all_cut_areas_match_geometry'))){throw 'Independent oblique incident power/fit hard checks failed'}
  $planes=@(Need $fit 'plane_data');$loo=@(Need $fit 'leave_one_out');if($planes.Count -ne 5 -or $loo.Count -ne 5){throw 'All five fixed independent electric fit planes and LOO required'}
  foreach($p in $planes){if(!(Is-True (Need $p 'dephased'))){throw 'Bloch transverse phase was not removed before projection'}}
  foreach($p in $loo){if([Math]::Abs((Num $p 'Pin_W')-1) -gt .005 -or (Num $p 'design_condition_2') -gt 10){throw 'Oblique source leave-one-out hard check failed'}}
  $basis=Need $Row 'expected_modes';$channels=@(Need $Metrics 'order_powers');$counts=@(Need $Metrics 'actual_independent_spatial_channel_counts')
  if($channels.Count -ne (Num $basis 'all_polarization_channels') -or $counts.Count -ne (Num $basis 'spatial_propagating_orders')){throw 'Physical mode inventory differs from independently enumerated lattice'}
  foreach($count in $counts){if($count -ne 2){throw 'A propagating spatial order does not have both polarizations'}}
}
function Assert-Plan {
 $script:plan=Json (Join-Path $stage 'plan.json');if($plan.stage_root -ne $stage -or (@($plan.cases.id) -join '|') -ne ($CaseIds -join '|')){throw 'Plan scope/stage changed'}
 $g=Need $plan 'geometry';foreach($p in @(@('pitch_um',80),@('opening_flat_to_flat_um',71.74019438588937),@('air_cap_um',10.5))){Equal-Number $g $p[0] $p[1]}
 foreach($id in $CaseIds){$row=Case-Row $id;$want=Case-Contract $id;foreach($n in @('height_um','opening_flat_to_flat_um','bottom_flat_to_flat_um','film_um','lambda_um','theta_deg','phi_deg')){Equal-Number $row $n $want.$n};if($row.polarization -ne $want.polarization -or $row.mesh -ne $want.mesh){throw 'Polarization/mesh changed'};Assert-IndependentModes $row
  $a=Analytic-Geometry $row;$e=Need $row 'expected_geometry';if($want.mesh -eq 'm4_air29_wft3_localtop24'){Equal-Number $e 'top_port_hmax_um' 2.4}elseif($null -ne (Need $e 'top_port_hmax_um')){throw 'Unapproved coarse-case top mesh refinement'};foreach($x in @(@('height_um',$want.height_um),@('top_plane_z_um',($want.height_um+10.5)),@('air_cap_um',10.5))){Equal-Number $e $x[0] $x[1]};if([Math]::Abs((Num $e 'film_volume_um3')/($a.film_volume_m3*1e18)-1) -gt 1e-10){throw 'Height-dependent film volume mismatch'}
 }
 $r=Need $plan 'resource_plan';foreach($p in @(@('cores_per_case',32),@('concurrency',1),@('owned_RSS_soft_guard_GiB',112),@('owned_RSS_hard_guard_GiB',128),@('Windows_free_before_admission_GiB_min',100),@('host_free_memory_runtime_minimum_GiB',64),@('host_CPU_before_admission_percent_max',75),@('actual_DOF_max',8000000),@('tetrahedra_max',2500000),@('min_tetrahedron_quality',.2),@('D_free_GiB_min',100))){Equal-Number $r $p[0] $p[1]}
 $q=Need $plan 'queue_plan';if($null -ne $q.new_admission_window_hours -or $q.time_window_enforced -ne $false){throw 'No time window authorized'};Equal-Number $q 'case_count_max' 2;if((@($q.fixed_order) -join '|') -ne ($CaseIds -join '|') -or (@($q.target_duration_hours) -join '|') -ne '6|16' -or $null -ne $q.per_case_wall_timeout_s -or $null -ne $q.cumulative_solver_budget_s){throw 'Finite window or no-timeout policy changed'}
 if(@($plan.mesh_only_case_ids).Count -ne 0 -or (@($plan.main_case_ids) -join '|') -ne ($MainCaseIds -join '|') -or (@($plan.conditional_case_ids) -join '|') -ne ($ConditionalCaseIds -join '|')){throw 'Closed main4 conditional2 scope changed'}
}
function Assert-Window {if((Hash $windowPath) -ne $script:windowHash){throw 'Context changed'};$w=Json $windowPath;if($w.purpose -ne 'G99_PAIR_UNWINDOWED_CONTEXT' -or $w.stage_root -ne $stage -or $w.manifest_sha256 -ne $manifestHash -or $w.time_window_enforced -ne $false -or $null -ne $w.new_admission_deadline_utc -or $w.started_utc -ne $startedUtc.ToString('o')){throw 'Unwindowed immutable context differs'}}
function Can-Admit([string]$Reason,[bool]$BeforePhysical=$false){Assert-Manifest;Assert-Window;if(Test-Path -LiteralPath (Join-Path $stage 'STOP_NEW_ADMISSION')){Save-Status 'STOPPED_NEW_ADMISSION_REQUESTED' $Reason;return $false};return $true}
function Assert-CalibrationBinding {
 $path=Flat-Path $script:CurrentCaseId;$a=Json $path;$phase=Source-Phase $script:CurrentCaseId;$ids=if($phase -eq 'Main'){$MainCaseIds}else{$ConditionalCaseIds}
 if($null -eq $script:calibrationReview -or $script:calibrationReview.path -ne $path -or (Hash $path) -ne $script:calibrationReview.sha256){throw 'Immutable phase-flat calibration receipt missing or changed'}
 if($a.state -ne (Flat-State $script:CurrentCaseId) -or $a.phase -ne $phase -or $a.manifest_sha256 -ne $manifestHash -or $a.queue_window_sha256 -ne $windowHash -or @($a.points).Count -ne $ids.Count){throw 'Native phase-source scope/window changed'}
 foreach($id in $ids){$r=Case-Contract $id;$points=@($a.points|Where-Object {$_.approved_case_id -eq $id});if($points.Count -ne 1){throw 'Required exactcase source proof absent'};$p=$points[0];if($p.state -ne 'PASS' -or $p.polarization -ne $r.polarization -or $p.wavelength_um -ne $r.lambda_um -or $p.theta_deg -ne $r.theta_deg -or $p.phi_global_deg -ne $r.phi_deg){throw 'Native phase-flat actual parameters differ'};foreach($h in $p.artifacts_sha256.PSObject.Properties){if([IO.Path]::GetFileName($h.Name) -ne $h.Name -or (Hash (Join-Path $p.directory $h.Name)) -ne $h.Value){throw 'Fresh phase-flat artifact changed'}}}
}
function Review-MeshRejection([string]$CaseId,[string]$OriginalError){
  Assert-Manifest;Assert-Window
  $path=Join-Path $stage ('controller_'+$CaseId+'_mesh.json')
  if(!(Test-Path -LiteralPath $path -PathType Leaf)){return $null}
  $raw=Json $path
  if($raw.state -ne 'STOPPED_REVIEW_REQUIRED' -or !(Has-Field $raw 'mesh_precheck_rejection') -or $null -eq $raw.mesh_precheck_rejection){return $null}
  $cs=Read-ControllerStatus $CaseId 'MeshOnly' 'STOPPED_REVIEW_REQUIRED';$r=Need $cs.status 'mesh_precheck_rejection'
  if((Need $r 'case_id') -ne $CaseId -or (Need $r 'mode') -ne 'MeshOnly' -or (Need $r 'kind') -notin @('DOF_CAP','MIN_QUALITY','TET_CAP') -or !(Is-False (Need $r 'solver_executed')) -or !(Is-False (Need $r 'automatic_coarsening_allowed')) -or !(Is-True (Need $r 'geometry_source_periodic_OOC_and_fingerprint_checks_passed'))){throw 'Typed mesh rejection not bound to safe unsolved preflight'}
  Equal-Number $r 'actual_DOF_max' 8000000;Equal-Number $r 'min_tetrahedron_quality_required' .2
  $expectedReason=if($r.kind -eq 'DOF_CAP'){'Actual mesh admission inventory exceeds bounds: degrees_of_freedom'}elseif($r.kind -eq 'TET_CAP'){'Actual mesh admission inventory exceeds bounds: tetrahedra'}else{'Actual tetrahedron quality below0.2'}
  if((Need $r 'reason') -ne $expectedReason -or (Need $cs.status 'reason') -ne $expectedReason -or $OriginalError -ne $expectedReason -or $null -ne $cs.status.active){throw 'Only exact mesh-precheck failure without cleanup/process errors is skippable'}
  $dir=Join-Path (Join-Path $stage $CaseId) 'mesh_preflight';$mesh=Json (Join-Path $dir 'mesh.json');$g=Json (Join-Path $dir 'geometry.json')
  $kind=Assert-ExactMesh $mesh $g $CaseId $true
  if($kind -ne $r.kind){throw 'Independent mesh-rejection classification differs'}
  foreach($n in @('degrees_of_freedom','tetrahedra','min_tetrahedron_quality')){if((Num $r $n) -ne (Num $mesh $n)){throw 'Typed rejection does not describe saved actual mesh'}}
  [void](Assert-Receipt $dir $CaseId 'MeshOnly')
  $exit=Json (Join-Path $dir 'native_exit_receipt.json');if((Need $exit 'case_id') -ne $CaseId -or (Num $exit 'compile_exit_code') -ne 0 -or (Num $exit 'native_exit_code') -ne 0){throw 'Only completed native mesh preflight can be skipped'}
  $hashes=[ordered]@{};foreach($n in @('mesh.json','geometry.json','case_receipt.json','linear_solver_override.json','mesh_only.mph','native_exit_receipt.json')){$hashes[$n]=Hash (Join-Path $dir $n)}
  $noJobs=Assert-NoStageNativeJobs;Assert-Manifest;Assert-Window
  return [ordered]@{case_id=$CaseId;film_um=(Case-Contract $CaseId).film_um;state='SKIPPED_MESH_RESOURCE_OR_QUALITY';accepted=$false;solver_executed=$false;new_physical_run=$false;rejection=$r;reason=$expectedReason;directory=$dir;controller_sha256=$cs.sha256;artifacts_sha256=$hashes;no_stage_native_jobs=$noJobs;no_retry_or_mesh_change=$true;not_completed_result=$true;recorded_utc=[DateTime]::UtcNow.ToString('o')}
}

function Mesh-RejectionKind($Dof,$Quality,$Tets){
 if(!(Number $Dof) -or !(Number $Quality) -or !(Number $Tets) -or $Dof -le 0 -or $Dof -ne [Math]::Floor($Dof) -or $Tets -le 0 -or $Tets -ne [Math]::Floor($Tets) -or $Quality -lt 0 -or $Quality -gt 1){throw 'Invalid mesh inventory is hard failure, not resource skip'}
 if($Dof -gt 8000000){return 'DOF_CAP'};if($Tets -gt 2500000){return 'TET_CAP'};if($Quality -lt .2){return 'MIN_QUALITY'};return $null
}
function Assert-OriginalMesh([double]$Thickness,$Tets,[string]$Sha){
 $row=Case-Row $script:CurrentCaseId;if($Thickness -ne $row.film_um -or $Sha -notmatch '^[a-f0-9]{64}$' -or !(Number $Tets) -or $Tets -le 0 -or $Tets -ne [Math]::Floor($Tets)){throw 'New case-bound actual mesh inventory missing'}
 $key=Mesh-FamilyKey $row;if($script:meshByThickness.ContainsKey($key)){$prior=$script:meshByThickness[$key];if($Sha -ne $prior.sha256 -or $Tets -ne $prior.tets){throw 'Matched s/p same-height same-profile actual mesh differs'}}
}
function Assert-ParentMeshProfile($Mesh,$Row){
 $short=(Need $Row 'mesh') -eq 'm4_air29_wft3_film02_short9';$variant=if($short){'short9_film02_Ti075_edge0225'}else{'original_m4_air29_wft3'}
 if((Need $Mesh 'mesh') -ne (Need $Row 'mesh') -or (Need $Mesh 'mesh_local_variant') -ne $variant -or !(Is-True (Need $Mesh 'mesh_local_top_boundary_size_readback_verified')) -or $null -ne (Need $Mesh 'mesh_hLocalTop_um') -or !(Is-True (Need $Mesh 'case_bound_mesh_size_readback_verified'))){throw 'Case-bound dynamic exact mesh-size readback missing'}
 $e=Need $Row 'expected_geometry';foreach($x in @(@('mesh_hAir_um','air_hmax_um'),@('mesh_hTi_um','Ti_hmax_um'),@('mesh_hEdge_um','edge_hmax_um'),@('mesh_hFilm_um','film_hmax_um'),@('film_hmax_um','film_hmax_um'))){Equal-Number $Mesh $x[0] (Num $e $x[1])}
 $p=Need $Mesh 'loss_region_mesh_admission';if((Need $p 'profile') -ne $Row.mesh -or !(Is-True (Need $p 'actual_profile_size_readback_verified')) -or !(Is-True (Need $p 'conditional_loss_mesh_admission')) -or (Need $p 'conditional9') -isnot [bool] -or [bool]$p.conditional9 -ne $short){throw 'Actual profile/loss-region mesh admission missing'}
 foreach($x in @(@('actual_air_hmax_um','air_hmax_um'),@('actual_Ti_boundary_hmax_um','Ti_hmax_um'),@('actual_edge_hmax_um','edge_hmax_um'),@('actual_film_hmax_um','film_hmax_um'))){Equal-Number $p $x[0] (Num $e $x[1])}
 $decay=(Num $Row 'lambda_um')/(4*[Math]::PI*(Num $Row 'k_film'));Equal-Number $p 'film_intensity_decay_length_um' $decay;Equal-Number $p 'film_hmax_over_intensity_decay' ((Num $e 'film_hmax_um')/$decay);Equal-Number $p 'conditional9_maximum_film_hmax_over_intensity_decay' .5
 if($short -and ((Num $Row 'lambda_um') -ne 9 -or (Num $p 'film_hmax_over_intensity_decay') -gt .5)){throw 'Conditional9 film attenuation-length mesh guard failed'}
}
function Assert-ExactMesh($Mesh,$Geometry,[string]$CaseId,[bool]$AllowResourceRejection=$false){
  $script:CurrentCaseId=$CaseId;$c=Case-Contract $CaseId;$row=Case-Row $CaseId;$thickness=[double]$c.film_um
  if((Need $Mesh 'case_id') -ne $CaseId -or (Need $Geometry 'case_id') -ne $CaseId -or (Need $Mesh 'state') -ne 'MESH_PREFLIGHT_ONLY' -or !(Is-False (Need $Mesh 'solver_executed')) -or (Need $Mesh 'mesh') -ne $row.mesh){throw 'Exact case-owned unsolved mesh required'}
  foreach($n in @('geometry_verified','periodic_mesh_verified','second_order_curl_elements_requested')){if(!(Is-True (Need $Mesh $n))){throw ('Mesh flag failed: '+$n)}}
  if(!(Is-True (Need $Geometry 'geometry_verified')) -or (Need $Geometry 'geometry') -ne $c.geometry_marker -or (Need $Geometry 'surface') -ne 'TI_IBC_WITH_SIO2_FLOOR_FILL' -or !(Is-True (Need $Geometry 'not_conformal_wall_or_land_coating')) -or !(Is-True (Need $Geometry 'single_point_spectral_angular_exploration')) -or (Need $Geometry 'proof_schema') -ne 'ANGLE_BAND_SOURCE_POWER_V2'){throw 'Fixed measured geometry/scope mismatch'}
  Assert-ProfileBinding $Geometry $row;Equal-Number $Geometry 'opening_flat_to_flat_um' $row.opening_flat_to_flat_um;Equal-Number $Geometry 'bottom_flat_to_flat_um' $row.bottom_flat_to_flat_um;Equal-Number $Mesh 'bottom_flat_to_flat_um' $row.bottom_flat_to_flat_um;$analytic=Analytic-Geometry $row;foreach($n in $analytic.Keys){if([Math]::Abs((Num $Geometry $n)/[double]$analytic[$n]-1) -gt 1e-6){throw ('Independent analytic geometry mismatch: '+$n)}}
  $expect=Need $row 'expected_geometry';$map=@{film_top_apothem_um=@('film_top_apothem_um',1.);film_volume_m3=@('film_volume_um3',1e-18);air_volume_m3=@('air_volume_um3',1e-18);total_computational_volume_m3=@('total_dielectric_volume_um3',1e-18);film_air_interface_area_m2=@('film_air_interface_area_um2',1e-12);top_area_m2=@('top_area_um2',1e-12)}
  foreach($n in $map.Keys){$bound=$map[$n];if([Math]::Abs((Num $Geometry $n)/((Num $expect $bound[0])*[double]$bound[1])-1) -gt 1e-6){throw 'Measured geometry differs from frozen per-case dimensions'}}
  if((Need $Mesh 'reported_dof_kind') -ne 'actual' -or [string](Need $Mesh 'electric_curl_element_order') -ne '2' -or (Need $Mesh 'mesh_fingerprint_scheme') -ne 'SHA256_RAW_COORDINATE_BITS_TET_TRI_CONNECTIVITY_ENTITY_IDS_V1' -or (Need $Mesh 'mesh_fingerprint_sha256') -notmatch '^[a-f0-9]{64}$'){throw 'Actual quadratic mesh fingerprint missing'}
  $dof=Num $Mesh 'degrees_of_freedom';$tets=Num $Mesh 'tetrahedra';$quality=Num $Mesh 'min_tetrahedron_quality'
  if($dof -le 0 -or $dof -ne [Math]::Floor($dof) -or $tets -le 0 -or $tets -ne [Math]::Floor($tets) -or $quality -lt 0 -or $quality -gt 1){throw 'Invalid actual mesh inventory or unchanged tetrahedron cap'}
  $resourceKind=Mesh-RejectionKind $dof $quality $tets
  $modes=Need $row 'expected_modes';foreach($x in @(@('propagating_spatial_orders',(Num $modes 'spatial_propagating_orders')),@('polarization_channels',(Num $modes 'propagating_polarization_channels')),@('total_port_channels_including_evanescent',(Num $modes 'all_polarization_channels')),@('film_hmax_um',(Num $row.expected_geometry 'film_hmax_um')),@('mesh_hAir_um',2.9),@('mesh_hTi_um',(Num $row.expected_geometry 'Ti_hmax_um')),@('mesh_hEdge_um',(Num $row.expected_geometry 'edge_hmax_um')),@('mesh_hFilm_um',(Num $row.expected_geometry 'film_hmax_um')))){Equal-Number $Mesh $x[0] $x[1]}
  Equal-Number $Mesh 'floor_fill_thickness_um' $thickness
  Assert-OriginalMesh $thickness $tets (Need $Mesh 'mesh_fingerprint_sha256')
  $key=Mesh-FamilyKey $row
  if($script:meshByThickness.ContainsKey($key)){$prior=$script:meshByThickness[$key];if((Need $Mesh 'mesh_fingerprint_sha256') -ne $prior.sha256 -or (Num $Mesh 'tetrahedra') -ne $prior.tets){throw 'Same-thickness s/p exact mesh mismatch'}}
  $pairs=@(Need $Mesh 'periodic_mesh_pairs');if($pairs.Count -ne 3){throw 'Three actual periodic mesh pairs required'}
  foreach($p in $pairs){if(!(Is-True (Need $p 'canonical_triangle_multiset_equal')) -or (Num $p 'source_triangles') -le 0 -or (Num $p 'source_triangles') -ne (Num $p 'destination_triangles')){throw 'Actual periodic mesh mismatch'}}
  foreach($pair in @(@('wavelength_um','lambda_um'),@('theta_deg','theta_deg'),@('phi_global_deg','phi_deg'))){Equal-Number $Geometry $pair[0] (Num $row $pair[1])}
  $e=@(Need $Geometry 'global_E_direction');$want=@((Controller-Incidence $row).E_in)
  if($e.Count -ne 3){throw 'Case-specific true p source vector missing'};for($j=0;$j -lt 3;$j++){if(!(Number $e[$j]) -or [Math]::Abs([double]$e[$j]-$want[$j]) -gt 1e-12){throw 'Case-specific true p source vector changed'}}
  Assert-Ooc (Need $Mesh 'linear_solver_override') $CaseId
  if($null -ne $resourceKind -and !$AllowResourceRejection){throw ('Mesh preflight requires case skip: '+$resourceKind)}
  foreach($o in @($Mesh,$Geometry)){Equal-Number $o 'height_um' $row.height_um;Equal-Number $o 'air_cap_um' 10.5;Equal-Number $o 'top_plane_z_um' ($row.height_um+10.5)}
  [void](Assert-PresolveProof (Need $Mesh 'presolve_port_wavevector_proof') $row);Assert-ParentMeshProfile $Mesh $row
  return $resourceKind
}
function Assert-Receipt([string]$Directory,[string]$CaseId,[string]$Mode){$r=Original-Assert-Receipt $Directory $CaseId $Mode;$p=Join-Path $Directory 'presolve_port_wavevectors.json';if((Need $r 'presolve_port_wavevector_proof_sha256') -ne (Hash $p)){throw 'Receipt port proof changed'};[void](Assert-PresolveProof (Json $p) (Case-Row $CaseId));return $r}
function Review-MeshAndAdmit([string]$CaseId) {
    $r=Original-Review-MeshAndAdmit $CaseId
    $mesh=Json (Join-Path $r.mesh_directory 'mesh.json')
    $key=Mesh-FamilyKey (Case-Row $CaseId)
    $script:meshByThickness[$key]=[ordered]@{case_id=$CaseId;sha256=$mesh.mesh_fingerprint_sha256;tets=$mesh.tetrahedra}
    return $r
}

function Review-PhysicalResult([string]$CaseId,$Admission,[DateTime]$CaseStarted){$r=Original-Review-PhysicalResult $CaseId $Admission $CaseStarted;$m=Json (Join-Path $r.directory 'metrics.json');Equal-Number $m 'height_um' (Case-Contract $CaseId).height_um;Equal-Number $m 'bottom_flat_to_flat_um' (Case-Contract $CaseId).bottom_flat_to_flat_um;$p=Json (Join-Path $r.directory 'presolve_port_wavevectors.json');if(($p|ConvertTo-Json -Depth 60 -Compress) -ne ($m.presolve_port_wavevector_proof|ConvertTo-Json -Depth 60 -Compress)){throw 'Physical metrics proof differs'};$r.artifacts_sha256['presolve_port_wavevectors.json']=Hash (Join-Path $r.directory 'presolve_port_wavevectors.json');$r['height_um']=(Case-Contract $CaseId).height_um;return $r}

function Save-Status([string]$State,[string]$Reason=$null) {
    if (!$script:statusOwned) { return }
    Write-Json $statusPath ([ordered]@{
        schema=1; purpose='G99_P80_MATCHED_WALL_FULLWAVE_PILOT'; state=$State; reason=$Reason
        stage_root=$stage; started_utc=$startedUtc.ToString('o'); updated_utc=[DateTime]::UtcNow.ToString('o')
        finished_utc=$finishedUtc; active=$active; resource=$resourceSnapshot
        planned_case_ids=$CaseIds; mesh_preflights=@($preflight.ToArray()); completed=@($completed.ToArray())
        skipped=@($skipped.ToArray()); calibration=$calibrationReview; pair_mesh_admission_sha256=$pairAdmissionHash
        pair_comparison=$pairComparison; peak_owned_rss_bytes=$peakRss
        manifest_sha256=$manifestHash; queue_window_sha256=$windowHash
        estimated_completion=[ordered]@{
            pair_solver_hours_after_admission=@(6,16); includes_resource_wait=$false; includes_model_or_mesh_time=$false
        }
        concurrency=1; cores_per_case=32; automatic_retries=0; solver_wall_timeout_s=$null
        formal_qualification=$false; mesh_convergence_established=$false; angular_directionality_established=$false
    })
}

function Write-Comparison {$p=Join-Path $stage 'exploration_comparison.json';Write-Json $p ([ordered]@{schema=1;case_rows=@($completed.ToArray());skipped=@($skipped.ToArray());saved_diagnostic=$diagnosticReview;shortwave_precheck=$shortwaveReview;warning_rows_are_not_accepted=$true;formal_qualification=$false;mesh_convergence_established=$false;complete_angular_directionality_established=$false;original_results_modified=$false});$script:pairComparison=[ordered]@{path=$p;sha256=Hash $p;rows=$completed.Count}}
function Test-LocalPolicy {
 foreach($id in $CaseIds){$r=Case-Contract $id;$inc=Controller-Incidence $r;$dot=0.;for($j=0;$j -lt 3;$j++){$dot+=$inc.k_in[$j]*$inc.E_in[$j]};if([Math]::Abs($dot) -gt 1e-8){throw 'Source not transverse'}}
 foreach($x in @(@(8000000,.2,2500000,$null),@(8000001,.2,100,'DOF_CAP'),@(5000000,.19,100,'MIN_QUALITY'),@(5000000,.2,2500001,'TET_CAP'))){if((Mesh-RejectionKind $x[0] $x[1] $x[2]) -ne $x[3]){throw 'Unchanged mesh guard fixture failed'}}
 return [ordered]@{state='PASS_G99_QUEUE_POLICY_SELF_TEST';physical_case_count=2;primary=2;conditional=0;cores_per_case=32;concurrency=1;time_window_enforced=$false;solver_executed=$false;remote_writes=$false}
}
function Start-PhaseFlat([string]$Phase){
 $script:CurrentCaseId=if($Phase -eq 'Main'){$MainCaseIds[0]}else{$ConditionalCaseIds[0]};Save-Status ('STARTING_'+$Phase.ToUpperInvariant()+'_NATIVE_FLAT_SOURCE_CONTROLS')
 $flatText=(& (Join-Path $stage 'run_flat_controls.ps1') -Phase $Phase -ExpectedManifestSha $manifestHash -ExpectedWindowSha $windowHash|Out-String).Trim();$script:calibrationReview=$flatText|ConvertFrom-Json
 if($calibrationReview.state -ne (Flat-State $script:CurrentCaseId)){throw 'Phase fresh flat native source controls failed'};Assert-CalibrationBinding
}




function Wait-QueueResources([string]$NextPhase) {
    while ($true) {
        Assert-Manifest
        Assert-Window
        if (Test-Path -LiteralPath (Join-Path $stage 'STOP_NEW_ADMISSION')) {
            Save-Status 'STOPPED_NEW_ADMISSION_REQUESTED' $NextPhase
            return $false
        }
        [void](Assert-NoStageNativeJobs)
        $os=Get-CimInstance Win32_OperatingSystem
        $cpu=(Get-CimInstance Win32_Processor | Measure-Object LoadPercentage -Average).Average
        $volume=Get-Volume -DriveLetter D
        $healthy=$volume.HealthStatus -eq 'Healthy' -and
            @($volume.OperationalStatus | Where-Object { [string]$_ -ne 'OK' }).Count -eq 0 -and
            @($volume.OperationalStatus).Count -gt 0
        if ($null -eq $os.FreePhysicalMemory -or $null -eq $cpu -or $null -eq $volume.SizeRemaining) {
            throw 'Resource counters unavailable for G99 admission'
        }
        $script:resourceSnapshot=[ordered]@{
            verified_utc=[DateTime]::UtcNow.ToString('o'); free_gib=[double]$os.FreePhysicalMemory/1MB
            cpu_percent=[double]$cpu; D_free_gib=[double]$volume.SizeRemaining/1GB
            D_healthy=$healthy; next_phase=$NextPhase
        }
        $reason=@()
        if ($resourceSnapshot.free_gib -lt 100) { $reason+='HOST_FREE_MEMORY_BELOW_100_GiB' }
        if ($resourceSnapshot.cpu_percent -ge 75) { $reason+='HOST_CPU_AT_OR_ABOVE_75_PERCENT' }
        if (!$resourceSnapshot.D_healthy -or $resourceSnapshot.D_free_gib -lt 100) { $reason+='D_VOLUME_UNHEALTHY_OR_BELOW_100_GiB' }
        if ($reason.Count -eq 0) { return $true }
        Save-Status 'WAITING_RESOURCES' ($NextPhase+': '+($reason -join ', '))
        Start-Sleep -Seconds 60
    }
}
function Test-Prerequisites {
    Assert-Manifest
    Assert-Plan
    foreach ($id in $CaseIds) {
        foreach ($name in @($id, ('admission_'+$id+'.json'), ('controller_'+$id+'_mesh.json'), ('controller_'+$id+'_solve.json'))) {
            if (Test-Path -LiteralPath (Join-Path $stage $name)) { throw 'G99 case previously attempted; preserve its artifacts' }
        }
    }
    [void](Assert-NoStageNativeJobs)
    return [ordered]@{
        state='PASS_G99_LIVE_IMMUTABLE_PREREQUISITES'; physical_case_count=2
        both_meshes_required_before_solves=$true; fresh_phase_flats_required=$true
        solver_executed=$false; remote_writes=$false
    }
}
if ($QueueSelfTest) { Test-LocalPolicy | ConvertTo-Json -Depth 20; return }
if ($PrerequisiteSelfTest) { Test-Prerequisites | ConvertTo-Json -Depth 30; return }

$lock=[IO.File]::Open((Join-Path $stage 'queue.lock'),[IO.FileMode]::OpenOrCreate,[IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
try {
    if ((Test-Path -LiteralPath $statusPath) -or (Test-Path -LiteralPath $windowPath)) { throw 'Existing G99 context; preserve previous attempt' }
    $pre=Test-Prerequisites
    $startedUtc=[DateTime]::UtcNow
    $window=[ordered]@{
        schema=1; purpose='G99_PAIR_UNWINDOWED_CONTEXT'; stage_root=$stage
        manifest_sha256=$manifestHash; started_utc=$startedUtc.ToString('o')
        new_admission_deadline_utc=$null; time_window_enforced=$false
        automatic_retries=0; concurrency=1; maximum_case_count=2; both_meshes_before_solves=$true
    }
    New-Json $windowPath $window
    $windowHash=Hash $windowPath
    $statusOwned=$true
    Save-Status 'STARTING_G99_PAIR'

    if (!(Wait-QueueResources 'native flat source controls')) { return }
    Start-PhaseFlat 'Main'
    Write-Comparison

    foreach ($id in $CaseIds) {
        $script:CurrentCaseId=$id
        $script:active=$null
        if (!(Wait-QueueResources ('mesh preflight '+$id))) { return }
        try {
            Invoke-Controller $id 'MeshOnly'
            $admission=Review-MeshAndAdmit $id
        } catch {
            $skip=Review-MeshRejection $id $_.Exception.Message
            if ($null -eq $skip) { throw }
            $skipped.Add($skip)
            $script:active=$null
            Save-Status 'MESH_PREFLIGHT_REJECTED' $skip.reason
            continue
        }
        $review=(Json $admission.path).review
        $entry=[ordered]@{
            case_id=$id; admission_sha256=$admission.sha256
            mesh_sha256=$admission.mesh_sha256; mesh_model_sha256=$admission.mesh_model_sha256
            actual_DOF=$review.actual_DOF; tetrahedra=$review.tetrahedra; minimum_quality=$review.minimum_quality
            mesh_fingerprint=$review.exact_mesh_fingerprint; solver_executed=$false
        }
        $preflight.Add($entry)
        $script:admissions[$id]=$admission
        $script:active=$null
        Save-Status 'MESH_PREFLIGHT_PASSED'
    }
    if ($preflight.Count -ne $CaseIds.Count) {
        $finishedUtc=[DateTime]::UtcNow.ToString('o')
        Write-Comparison
        Save-Status 'STOPPED_PAIR_MESH_ADMISSION' 'Both G99 meshes must satisfy the frozen preflight before either physical solve'
        return
    }
    $pairPath=Join-Path $stage 'pair_mesh_admission.json'
    New-Json $pairPath ([ordered]@{
        state='PASS_G99_TWO_MESH_PREFLIGHTS'; stage_root=$stage
        manifest_sha256=$manifestHash; queue_window_sha256=$windowHash; entries=@($preflight.ToArray())
        reviewed_utc=[DateTime]::UtcNow.ToString('o'); physical_solver_executed=$false
    })
    $script:pairAdmissionHash=Hash $pairPath
    Save-Status 'BOTH_MESH_PREFLIGHTS_PASSED'

    foreach ($id in $CaseIds) {
        $script:CurrentCaseId=$id
        if (!(Wait-QueueResources ('physical solve '+$id))) { return }
        if ((Hash $pairPath) -ne $pairAdmissionHash) { throw 'Two-mesh preflight receipt changed' }
        $caseStarted=[DateTime]::UtcNow
        Invoke-Controller $id 'SolveApproved'
        $result=Review-PhysicalResult $id $admissions[$id] $caseStarted
        $completed.Add($result)
        $script:active=$null
        Write-Comparison
        Save-Status 'CASE_PHYSICS_COMPLETE'
    }
    $finishedUtc=[DateTime]::UtcNow.ToString('o')
    Write-Comparison
    Save-Status 'COMPLETE_G99_NORMAL_PAIR'
} catch {
    $reason=$_.Exception.Message
    $script:active=$null
    if ($statusOwned) {
        $finishedUtc=[DateTime]::UtcNow.ToString('o')
        try { Write-Comparison; Save-Status 'STOPPED_REVIEW_REQUIRED' $reason } catch {}
    }
    throw
} finally {
    if ($null -ne $lock) { $lock.Dispose() }
}
