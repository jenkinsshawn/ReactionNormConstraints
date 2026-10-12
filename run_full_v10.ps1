<#
.SYNOPSIS
  Launches the full ReactionNormConstraints v10.0.0 reruns with Rscript.

.DESCRIPTION
  Runs run_full_v10.R from the folder that holds this script (the repository
  root). Output goes to .\sim_output_v10\ and .\audit_output_v10\; per-part logs
  and a run summary go to .\logs\. Expect about two to three hours at full
  scale; runtime estimates are printed at start-up.

.PARAMETER Parts
  Comma-separated subset of: sim_v9, sim_v10, audit, check. Default: all four.

.PARAMETER Smoke
  Run in SMOKE_TEST mode (replicate counts scaled to ~2%) as a pipeline check.

.PARAMETER Cores
  Number of parallel R workers. Default: 85% of logical processors.

.PARAMETER Rscript
  Full path to Rscript.exe if it is not on PATH.

.PARAMETER DryRun
  Print the plan and runtime estimates without running anything.

.EXAMPLE
  .\run_full_v10.ps1 -Smoke
.EXAMPLE
  .\run_full_v10.ps1 -Parts "sim_v10,audit" -Cores 24
.NOTES
  If script execution is blocked, run once in this window:
    Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
#>
param(
  [string]$Parts = "sim_v9,sim_v10,audit,check",
  [switch]$Smoke,
  [int]$Cores = 0,
  [string]$Rscript = "",
  [switch]$DryRun
)
$ErrorActionPreference = "Stop"
if ((Get-Location).ProviderPath -ne $PSScriptRoot) { Set-Location -LiteralPath $PSScriptRoot }

if (-not $Rscript) {
  $cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue
  if ($cmd) {
    $Rscript = $cmd.Source
  } else {
    # Registry entries written by the CRAN installer, then the default install folders
    foreach ($key in "HKCU:\SOFTWARE\R-core\R", "HKLM:\SOFTWARE\R-core\R") {
      $ip = (Get-ItemProperty -Path $key -Name InstallPath -ErrorAction SilentlyContinue).InstallPath
      if ($ip -and (Test-Path -LiteralPath (Join-Path $ip "bin\Rscript.exe"))) { $Rscript = Join-Path $ip "bin\Rscript.exe"; break }
    }
    if (-not $Rscript) {
      $cand = Get-ChildItem -Path "$env:ProgramFiles\R\R-*\bin\Rscript.exe", "$env:LOCALAPPDATA\Programs\R\R-*\bin\Rscript.exe", "$env:USERPROFILE\AppData\Local\Programs\R\R-*\bin\Rscript.exe" -ErrorAction SilentlyContinue |
        Sort-Object -Property FullName -Descending | Select-Object -First 1
      if ($cand) { $Rscript = $cand.FullName }
    }
  }
}
if (-not $Rscript -or -not (Test-Path -LiteralPath $Rscript)) {
  throw "Rscript.exe not found. Pass it explicitly, for example: -Rscript 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'"
}

$rargs = @("run_full_v10.R", "--parts=$Parts")
if ($Smoke)   { $rargs += "--smoke" }
if ($Cores -gt 0) { $rargs += "--cores=$Cores" }
if ($DryRun)  { $rargs += "--dry-run" }

Write-Host "Rscript : $Rscript"
Write-Host "Folder  : $PSScriptRoot"
Write-Host "Command : Rscript $($rargs -join ' ')"
$sw = [System.Diagnostics.Stopwatch]::StartNew()
# R writes warnings to stderr; do not let PowerShell treat them as terminating errors
$ErrorActionPreference = "Continue"
& $Rscript @rargs
$code = $LASTEXITCODE
$sw.Stop()
Write-Host ("Finished with exit code {0} after {1:N1} minutes." -f $code, $sw.Elapsed.TotalMinutes)
exit $code
