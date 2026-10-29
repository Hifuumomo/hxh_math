$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    for ($pass = 1; $pass -le 3; $pass++) {
        & xelatex -interaction=nonstopmode -halt-on-error main.tex
        if ($LASTEXITCODE -ne 0) {
            throw "XeLaTeX 编译失败，请检查 main.log。"
        }
    }
} finally {
    Pop-Location
}
