@echo off
setlocal
cd /d "%~dp0"
python -m kl_convergence.noise_sensitivity ^
  --config configs\server_kl_convergence_v1_0_1.yaml ^
  --frozen-dir ..\..\..\data\inputs ^
  --baseline-dir ..\..\..\data\kl_convergence_v1_0_1 ^
  --out results_noise_sensitivity_primary_dt_v1_0_1
if errorlevel 1 exit /b %errorlevel%
endlocal
