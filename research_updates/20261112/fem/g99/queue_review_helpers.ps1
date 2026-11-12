# Extracted from the pinned queue; only the G99 evidence bindings change.
function Has-Field($Object,[string]$Name){return $null -ne $Object -and $null -ne $Object.PSObject.Properties[$Name]}

function Need($Object,[string]$Name){if(!(Has-Field $Object $Name)){throw ('Missing queue evidence field: '+$Name)};return $Object.$Name}

function Is-True($Value){return $Value -is [bool] -and $Value}

function Is-False($Value){return $Value -is [bool] -and !$Value}

function Number($Value){try{$n=[double]$Value}catch{return $false};return $null -ne $Value -and $Value -isnot [bool] -and ![double]::IsNaN($n) -and ![double]::IsInfinity($n)}

function Num($Object,[string]$Name){$v=Need $Object $Name;if(!(Number $v)){throw ('Nonfinite/non-numeric queue evidence: '+$Name)};return [double]$v}

function Equal-Number($Object,[string]$Name,[double]$Expected){if([Math]::Abs((Num $Object $Name)-$Expected) -gt 1e-12){throw ('Fixed queue evidence mismatch: '+$Name)}}

function Hash([string]$Path){if(!(Test-Path -LiteralPath $Path -PathType Leaf) -or (Get-Item -LiteralPath $Path).Length -le 0){throw ('Missing/empty queue artifact: '+$Path)};return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()}

function Json([string]$Path){return [IO.File]::ReadAllText($Path,[Text.Encoding]::UTF8)|ConvertFrom-Json}

function Write-Json([string]$Path,$Value){$tmp=$Path+'.tmp';[IO.File]::WriteAllText($tmp,($Value|ConvertTo-Json -Depth 55),[Text.UTF8Encoding]::new($false));Move-Item -LiteralPath $tmp -Destination $Path -Force}

function New-Json([string]$Path,$Value){
  # CreateNew is the non-overwrite boundary, including admissions and window.
  $stream=[IO.File]::Open($Path,[IO.FileMode]::CreateNew,[IO.FileAccess]::Write,[IO.FileShare]::Read)
  try{$bytes=[Text.UTF8Encoding]::new($false).GetBytes(($Value|ConvertTo-Json -Depth 55));$stream.Write($bytes,0,$bytes.Length);$stream.Flush()}finally{$stream.Dispose()}
}

function Beijing([DateTime]$Utc){return $Utc.ToUniversalTime().AddHours(8).ToString("yyyy-MM-dd HH:mm:ss 'UTC+08:00'",[Globalization.CultureInfo]::InvariantCulture)}

function Case-Row([string]$Id){$rows=@($script:plan.cases|Where-Object {$_.id -eq $Id});if($Id -notin $CaseIds -or $rows.Count -ne 1){throw 'Frozen case binding missing'};return $rows[0]}

function Original-Assert-Receipt([string]$Directory,[string]$CaseId,[string]$Mode){
  $r=Json (Join-Path $Directory 'case_receipt.json')
  if((Need $r 'case_id') -ne $CaseId -or (Need $r 'mode') -ne $Mode -or (Need $r 'attempt_id') -ne 'initial' -or (Need $r 'manifest_sha256') -ne $script:manifestHash -or (Need $r 'queue_window_sha256') -ne $script:windowHash){throw 'Case receipt not bound to this unique queue invocation'}
  foreach($b in @(@('core_template_sha256','CavityFloorFilm.java'),@('controller_sha256','run_floorfilm_controller.ps1'),@('plan_sha256','plan.json'),@('material_sha256','material.json'))){if((Need $r $b[0]) -ne (Hash (Join-Path $stage $b[1]))){throw 'Receipt frozen source hash mismatch'}}
  foreach($b in @(@('generated_source_sha256','CaseEntry.java'),@('compiled_case_sha256','CaseEntry.class'))){if((Need $r $b[0]) -ne (Hash (Join-Path $Directory $b[1]))){throw 'Receipt generated/compiled case hash mismatch'}}
  $oocPath=Join-Path $Directory 'linear_solver_override.json'
  if((Need $r 'linear_solver_override_sha256') -ne (Hash $oocPath)){throw 'Receipt OOC readback artifact changed'}
  Assert-Ooc (Json $oocPath) $CaseId
  if($Mode -eq 'MeshOnly' -and !(Is-True (Need $r 'mesh_only'))){throw 'Mesh receipt is not mesh-only'}
  if($Mode -eq 'SolveApproved' -and !(Is-False (Need $r 'mesh_only'))){throw 'Physical receipt is still mesh-only'}
  return $r
}

function Read-ControllerStatus([string]$CaseId,[string]$Mode,[string]$ExpectedState){
  $short=if($Mode -eq 'MeshOnly'){'mesh'}else{'solve'}
  $path=Join-Path $stage ('controller_'+$CaseId+'_'+$short+'.json');$s=Json $path
  if(($ExpectedState -split '\|') -notcontains [string](Need $s 'state') -or (Need $s 'case_id') -ne $CaseId -or (Need $s 'mode') -ne $Mode -or (Need $s 'attempt_id') -ne 'initial' -or (Need $s 'manifest_sha256') -ne $script:manifestHash -or (Need $s 'queue_window_sha256') -ne $script:windowHash){throw ('Controller did not reach required bound terminal state: '+$ExpectedState)}
  $rss=Num $s 'peak_owned_rss_bytes';if($rss -lt 0){throw 'Invalid controller peak RSS'};if($rss -gt $script:peakRss){$script:peakRss=[long]$rss}
  return [pscustomobject]@{path=$path;status=$s;sha256=Hash $path}
}

function Invoke-Controller([string]$CaseId,[string]$Mode){
  Assert-Manifest;Assert-Window;Assert-CalibrationBinding
  $script:active=[ordered]@{case_id=$CaseId;phase=$Mode;started_utc=[DateTime]::UtcNow.ToString('o');controller_status_path=Join-Path $stage ('controller_'+$CaseId+'_'+$(if($Mode -eq 'MeshOnly'){'mesh'}else{'solve'})+'.json')}
  Save-Status ('RUNNING_'+$Mode.ToUpperInvariant())
  # Synchronous invocation is intentional: the child owns resource protection.
  # The queue never sends a kill solely because its new-admission window ends.
  & $controller -Mode $Mode -CaseId $CaseId -ExpectedManifestSha $script:manifestHash -StageRoot $stage -ExpectedWindowSha $script:windowHash
  Assert-Manifest;Assert-Window
}

function Original-Review-MeshAndAdmit([string]$CaseId){
  $cs=Read-ControllerStatus $CaseId 'MeshOnly' 'COMPLETE_MESH_PREFLIGHT_AWAITING_EXPLICIT_ADMISSION'
  $dir=Join-Path (Join-Path $stage $CaseId) 'mesh_preflight'
  $c=Need $cs.status 'completed'
  if((Need $c 'directory') -ne $dir -or (Num $c 'exit_code') -ne 0 -or !(Is-True (Need $c 'eligible_for_manual_review')) -or !(Is-False (Need $c 'solver_executed'))){throw 'Completed mesh controller evidence inconsistent'}
  $mesh=Json (Join-Path $dir 'mesh.json');$geometry=Json (Join-Path $dir 'geometry.json')
  [void](Assert-ExactMesh $mesh $geometry $CaseId);[void](Assert-Receipt $dir $CaseId 'MeshOnly')
  $modelHash=Hash (Join-Path $dir 'mesh_only.mph')
  if((Need $c 'mesh_model_sha256') -ne $modelHash){throw 'Persisted mesh model changed after controller validation'}
  $entry=[ordered]@{case_id=$CaseId;allow_solve=$true;mesh_relative_directory=($CaseId+'\mesh_preflight');mesh_sha256=Hash (Join-Path $dir 'mesh.json');geometry_sha256=Hash (Join-Path $dir 'geometry.json');receipt_sha256=Hash (Join-Path $dir 'case_receipt.json');mesh_model_sha256=$modelHash;core_sha256=Hash (Join-Path $stage 'CavityFloorFilm.java');plan_sha256=Hash (Join-Path $stage 'plan.json');controller_sha256=Hash $controller;manifest_sha256=$script:manifestHash;presolve_port_wavevector_proof_sha256=Hash (Join-Path $dir 'presolve_port_wavevectors.json');fresh_flat_admission_sha256=$calibrationReview.sha256;fresh_flat_admission_path=Flat-Path $CaseId}
  $admission=[ordered]@{schema=1;purpose='FLOORFILM_CAVITY_MESH_ADMISSION';stage_root=$stage;queue_window_sha256=$script:windowHash;reviewed_by='Codex user-approved frozen finite exploration exact per-case mesh and OOC evidence review';reviewed_utc=[DateTime]::UtcNow.ToString('o');entries=@($entry);review=[ordered]@{actual_DOF=Num $mesh 'degrees_of_freedom';tetrahedra=Num $mesh 'tetrahedra';minimum_quality=Num $mesh 'min_tetrahedron_quality';exact_mesh_fingerprint=Need $mesh 'mesh_fingerprint_sha256';source_controller_assertions_passed=$true;queue_independent_exact_mesh_and_OOC_readback_passed=$true;only_floorfill_not_conformal=$true;solver_started=$false;source_definition='true s in global xz incidence plane';wavelength_um=(Case-Contract $CaseId).lambda_um;theta_deg=(Case-Contract $CaseId).theta_deg;phi_global_deg=0;flat_calibration_sha256=$script:calibrationReview.sha256;floor_fill_thickness_um=(Case-Contract $CaseId).film_um;physical_thresholds_unchanged=$true}}
  Assert-Manifest;Assert-Window
  $path=Join-Path $stage ('admission_'+$CaseId+'.json');New-Json $path $admission
  return [pscustomobject]@{path=$path;sha256=Hash $path;mesh_controller_sha256=$cs.sha256;mesh_directory=$dir;mesh_sha256=$entry.mesh_sha256;geometry_sha256=$entry.geometry_sha256;receipt_sha256=$entry.receipt_sha256;mesh_model_sha256=$modelHash;mesh_peak_owned_rss_bytes=Num $cs.status 'peak_owned_rss_bytes';mesh_elapsed_s=([DateTime]::UtcNow-[DateTimeOffset]::Parse($cs.status.started_utc,[Globalization.CultureInfo]::InvariantCulture).UtcDateTime).TotalSeconds}
}

function Assert-Disposition($Evaluation,$Metrics){
  $checks=Need $Evaluation 'checks';$properties=@($checks.PSObject.Properties)
  if($properties.Count -lt 40){throw 'Independent evaluation check inventory missing'}
  $failed=@($properties|Where-Object {!(Is-True $_.Value)}|ForEach-Object {$_.Name});$hard=@($failed|Where-Object {$_ -notin $AllowedWarnings});$warnings=@($failed|Where-Object {$_ -in $AllowedWarnings})
  if((@($Evaluation.failed_checks|Sort-Object) -join '|') -ne (@($failed|Sort-Object) -join '|') -or (@($Evaluation.hard_failed_checks|Sort-Object) -join '|') -ne (@($hard|Sort-Object) -join '|') -or (@($Evaluation.warning_checks|Sort-Object) -join '|') -ne (@($warnings|Sort-Object) -join '|')){throw 'Recorded failure/warning/hard lists differ from full checks'}
  if($hard.Count -ne 0 -or !(Is-True (Need $Evaluation 'exploration_continue_allowed'))){throw ('Hard physical failure stops queue: '+($hard -join ','))}
  $clean=$failed.Count -eq 0
  if((Need $Evaluation 'accepted') -isnot [bool] -or [bool]$Evaluation.accepted -ne $clean){throw 'Acceptance was relabelled instead of retaining failures'}
  $expected=if($clean){'COMPLETE_EXPLORATORY_POINT_ACCEPTED'}else{'COMPLETE_WITH_FLUX_WARNING'}
  if((Need $Evaluation 'disposition') -ne $expected){throw 'Disposition differs from full check classification'}
  foreach($n in @('formal_qualification','mesh_convergence_established','raw_metrics_changed')){if(!(Is-False (Need $Evaluation $n))){throw 'No raw metric changes, convergence or qualification allowed'}}
  if($null -ne $Metrics){
    $native=Need $Metrics 'checks';$nativeFailed=0
    foreach($p in $native.PSObject.Properties){$name='native_'+$p.Name;$flag=Need $checks $name;if($p.Value -isnot [bool] -or $flag -isnot [bool] -or [bool]$p.Value -ne [bool]$flag){throw 'Original Java checks altered in independent evaluation'};if(!$p.Value){$nativeFailed++;if($name -notin $AllowedWarnings){throw ('Native hard failure: '+$name)}}}
    $nativeState=if($nativeFailed -eq 0){'PASSED_FLOOR_FILL_CAVITY_DIAGNOSTIC_ONLY'}else{'FAILED_FLOOR_FILL_CAVITY_DIAGNOSTIC'};if((Need $Metrics 'state') -ne $nativeState){throw 'Original Java state relabelled'}
    $port=1-(Num $Metrics 'R_total');$loss=Num $Metrics 'A_loss_total';$ti=Num $Metrics 'A_Ti_H';$film=Num $Metrics 'A_film';$tiE=Num $Metrics 'A_Ti_E';$qh=Num $Metrics 'A_film_builtin_Qh'
    if([Math]::Abs($loss-$ti-$film) -gt 1e-10 -or [Math]::Abs($port-$loss) -gt .005 -or [Math]::Abs($port-$loss)/[Math]::Max([Math]::Abs($port),1e-12) -gt .02 -or [Math]::Abs($ti-$tiE) -gt .005 -or [Math]::Abs($ti-$tiE)/[Math]::Max([Math]::Abs($ti),1e-12) -gt .02 -or [Math]::Abs($film-$qh) -gt .005 -or [Math]::Abs($film-$qh)/[Math]::Max([Math]::Abs($film),1e-12) -gt .02 -or [Math]::Abs((Num $Metrics 'Pin_independent_internal_E_fit_W')-1) -gt .005){throw 'Independent parent hard loss/material/source accounting failed'}
  }
  return [pscustomobject]@{accepted=$clean;warning_checks=$warnings;hard_failed_checks=$hard;exploration_continue_allowed=$true;disposition=$expected}
}

function Original-Review-PhysicalResult([string]$CaseId,$Admission,[DateTime]$CaseStarted){
  $cs=Read-ControllerStatus $CaseId 'SolveApproved' 'COMPLETE_FLOOR_FILL_CAVITY_DIAGNOSTIC_ONLY|COMPLETE_FLOOR_FILL_WITH_FLUX_WARNING'
  $dir=Join-Path (Join-Path $stage $CaseId) 'solve';$c=Need $cs.status 'completed'
  if((Need $c 'directory') -ne $dir -or (Num $c 'exit_code') -ne 0 -or !(Is-True (Need $c 'solver_executed')) -or (Need $cs.status 'admission_sha256') -ne $Admission.sha256 -or (Hash $Admission.path) -ne $Admission.sha256){throw 'Completed solve/admission binding inconsistent'}
  $evaluation=Json (Join-Path $dir 'evaluation.json');$recorded=Need $c 'evaluation'
  if((Need $evaluation 'case_id') -ne $CaseId -or (Need $recorded 'case_id') -ne $CaseId -or ($evaluation|ConvertTo-Json -Depth 50 -Compress) -ne ($recorded|ConvertTo-Json -Depth 50 -Compress)){throw 'Persisted/recorded independent evaluation differs'}
  $mesh=Json (Join-Path $dir 'mesh.json');$geometry=Json (Join-Path $dir 'geometry.json');[void](Assert-ExactMesh $mesh $geometry $CaseId);[void](Assert-Receipt $dir $CaseId 'SolveApproved')
  $m=Json (Join-Path $dir 'metrics.json')
  if((Need $m 'case_id') -ne $CaseId -or !(Is-True (Need $m 'solver_executed')) -or !(Is-True (Need $m 'solver_completed')) -or !(Is-True (Need $m 'linear_solver_override_reverified_before_physics'))){throw 'Actual physical solve/OOC evidence missing'}
  $disposition=Assert-Disposition $evaluation $m;Assert-Ooc (Need $m 'linear_solver_override') $CaseId
  $expectedTerminal=if($disposition.accepted){'COMPLETE_FLOOR_FILL_CAVITY_DIAGNOSTIC_ONLY'}else{'COMPLETE_FLOOR_FILL_WITH_FLUX_WARNING'};if($cs.status.state -ne $expectedTerminal){throw 'Controller terminal state disagrees with preserved checks'}
  $want=Case-Contract $CaseId
  if((Need $m 'polarization') -ne $want.polarization -or (Num $m 'floor_fill_thickness_um') -ne $want.film_um -or (Num $m 'wavelength_um') -ne $want.lambda_um -or (Num $m 'theta_deg') -ne $want.theta_deg -or (Num $m 'phi_global_deg') -ne 0){throw 'Physical point/thickness/polarization changed'}
  Assert-ActualSource $m (Case-Row $CaseId)
  if((Need $m 'mesh_fingerprint_sha256') -ne $mesh.mesh_fingerprint_sha256 -or (Need $m 'reviewed_mesh_fingerprint_sha256') -ne $mesh.mesh_fingerprint_sha256 -or (Num $m 'degrees_of_freedom') -ne $mesh.degrees_of_freedom -or (Num $m 'tetrahedra') -ne $mesh.tetrahedra){throw 'Physical metrics not bound to admitted exact mesh'}
  $hashes=[ordered]@{};foreach($n in @('metrics.json','evaluation.json','solved_unevaluated.mph','model.mph','mesh_only.mph','mesh.json','geometry.json','case_receipt.json','linear_solver_override.json','solved_artifact.json','done.json','native_exit_receipt.json')){$hashes[$n]=Hash (Join-Path $dir $n)}
  if((Need $c 'metrics_sha256') -ne $hashes['metrics.json'] -or (Need $c 'model_sha256') -ne $hashes['model.mph']){throw 'Completed physical artifact hashes changed'}
  $metrics=[ordered]@{};foreach($n in @('R_total','A_port','A_loss_total','A_film','A_Ti_H','R_total_net_flux','A_flux','A_surface_Ti_H','A_Ti_E','flux_total_loss_absolute_mismatch','flux_total_loss_relative_mismatch','side_pair_power_mismatch_over_Pin','edge_associated_absorption')){$metrics[$n]=Num $m $n}
  return [ordered]@{case_id=$CaseId;film_um=$want.film_um;polarization=$want.polarization;wavelength_um=$want.lambda_um;theta_deg=$want.theta_deg;phi_global_deg=0;state=$(if($disposition.accepted){'EXPLORATORY_SINGLE_POINT_ACCEPTED_NOT_MESH_CONVERGED'}else{'COMPLETE_WITH_FLUX_WARNING_NOT_ACCEPTED'});accepted=$disposition.accepted;exploration_continue_allowed=$true;warning_checks=$disposition.warning_checks;hard_failed_checks=$disposition.hard_failed_checks;directory=$dir;finished_utc=[DateTime]::UtcNow.ToString('o');case_wall_elapsed_s=([DateTime]::UtcNow-$CaseStarted).TotalSeconds;solver_elapsed_s=Num $m 'solve_seconds';peak_owned_rss_bytes=[Math]::Max([double]$Admission.mesh_peak_owned_rss_bytes,(Num $cs.status 'peak_owned_rss_bytes'));mesh_peak_owned_rss_bytes=$Admission.mesh_peak_owned_rss_bytes;solve_peak_owned_rss_bytes=Num $cs.status 'peak_owned_rss_bytes';mesh_fingerprint_sha256=$mesh.mesh_fingerprint_sha256;admission_sha256=$Admission.sha256;mesh_controller_sha256=$Admission.mesh_controller_sha256;solve_controller_sha256=$cs.sha256;queue_window_sha256=$script:windowHash;manifest_sha256=$script:manifestHash;artifacts_sha256=$hashes;metrics=$metrics;evaluation=$evaluation;formal_qualification=$false;mesh_convergence_established=$false;angular_directionality_established=$false}
}

function Assert-NoStageNativeJobs {
  $all=@(Get-CimInstance Win32_Process -ErrorAction Stop);$owned=@{};$found=$true
  foreach($p in $all){if(([string]$p.CommandLine).IndexOf($stage,[StringComparison]::OrdinalIgnoreCase) -ge 0 -or ([string]$p.ExecutablePath).IndexOf($stage,[StringComparison]::OrdinalIgnoreCase) -ge 0){$owned[[int]$p.ProcessId]=$true}}
  while($found){$found=$false;foreach($p in $all){if($owned.ContainsKey([int]$p.ParentProcessId) -and !$owned.ContainsKey([int]$p.ProcessId)){$owned[[int]$p.ProcessId]=$true;$found=$true}}}
  $jobs=@($all|Where-Object {$owned.ContainsKey([int]$_.ProcessId) -and $_.Name -match '^(?:comsol[^.]*|javaw?|mpiexec(?:\.hydra)?|hydra_pmi_proxy)\.exe$'})
  if($jobs.Count){throw ('Mesh skip forbidden while this stage still owns native jobs: '+(@($jobs|ForEach-Object {$_.ProcessId}) -join ','))}
  return [ordered]@{checked_utc=[DateTime]::UtcNow.ToString('o');stage_root=$stage;native_job_count=0;query='Win32_Process stage path and descendant closure, native COMSOL Java MPI identities';no_process_killed=$true}
}
