$ErrorActionPreference = 'Stop'
$here = Join-Path (Split-Path -Parent $MyInvocation.MyCommand.Path) 'ldpcQaun'
Push-Location $here
try {
    & .\build.bat
    if ($LASTEXITCODE -ne 0) { throw "Build failed with exit code $LASTEXITCODE" }

    foreach ($code in @('18-original', '18-ds3', '36', '36-4-6')) {
        & .\DS_CSS_decoder_parallel.exe --code $code --self-test
        if ($LASTEXITCODE -ne 0) { throw "Self-test failed for $code" }
    }

    foreach ($code in @('18-original', '18-ds3', '36', '36-4-6')) {
        & .\DS_CSS_decoder_parallel.exe --code $code --p-data 0.01 --p-meas 0.01 --iterations 30 --order 2 --shots 20 --seed 123
        if ($LASTEXITCODE -ne 0) { throw "Monte Carlo run failed for $code" }
    }
}
finally { Pop-Location }
