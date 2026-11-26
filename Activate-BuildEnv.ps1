<#
.SYNOPSIS
Activate the x64 MSVC, CUDA, CMake/Ninja and vcpkg build environment.
.DESCRIPTION
Run from PowerShell: & "$([Environment]::GetFolderPath('Desktop'))\Activate-BuildEnv.ps1"
Or double-click the Build Environment desktop shortcut to open an activated shell.
Changes apply to this process and its children. Existing CMake caches retain their settings.
#>
[CmdletBinding()]
param(
    [string]$VcpkgRoot = 'C:\opt\vcpkg',
    [string]$CudaRoot = $env:CUDA_PATH
)

& {
    $ErrorActionPreference = 'Stop'
    $vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
    if (-not (Test-Path -LiteralPath $vswhere)) { throw 'Visual Studio Installer / vswhere.exe was not found.' }
    $vsRoot = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $vsRoot) { throw 'No Visual Studio installation with the MSVC x64 tools was found.' }
    if (-not (Test-Path -LiteralPath "$VcpkgRoot\scripts\buildsystems\vcpkg.cmake")) { throw "vcpkg toolchain not found under $VcpkgRoot. Use -VcpkgRoot to specify its location." }
    if (-not $CudaRoot) {
        $cudaBase = "$env:ProgramFiles\NVIDIA GPU Computing Toolkit\CUDA"
        $latestCuda = Get-ChildItem -LiteralPath $cudaBase -Directory -ErrorAction SilentlyContinue |
            Where-Object { $_.Name -match '^v\d+\.\d+$' } |
            Sort-Object { [version]$_.Name.Substring(1) } -Descending | Select-Object -First 1
        if ($latestCuda) { $CudaRoot = $latestCuda.FullName }
    }
    if (-not $CudaRoot -or -not (Test-Path -LiteralPath "$CudaRoot\bin\nvcc.exe")) { throw 'CUDA nvcc.exe was not found. Use -CudaRoot to specify the toolkit directory.' }

    & "$vsRoot\Common7\Tools\Launch-VsDevShell.ps1" -Arch amd64 -HostArch amd64 -SkipAutomaticLocation
    if (-not (Get-Command cl.exe -ErrorAction SilentlyContinue)) { throw 'MSVC environment initialization failed.' }

    $cmakeBase = "$vsRoot\Common7\IDE\CommonExtensions\Microsoft\CMake"
    $paths = @("$CudaRoot\bin", $VcpkgRoot, "$cmakeBase\CMake\bin", "$cmakeBase\Ninja") |
        Where-Object { Test-Path -LiteralPath $_ -PathType Container }
    $env:PATH = (($paths + @($env:PATH -split ';') | Where-Object { $_ } | Select-Object -Unique) -join ';')
    $env:VCPKG_ROOT = $VcpkgRoot
    $env:VCPKG_DEFAULT_TRIPLET = 'x64-windows'
    $env:VCPKG_DEFAULT_HOST_TRIPLET = 'x64-windows'
    $env:CMAKE_TOOLCHAIN_FILE = "$VcpkgRoot\scripts\buildsystems\vcpkg.cmake"
    $env:CMAKE_GENERATOR = 'Ninja'
    $env:CUDA_PATH = $CudaRoot
    $env:CUDACXX = "$CudaRoot\bin\nvcc.exe"
    $env:CUDAHOSTCXX = (Get-Command cl.exe).Source
    $env:CC = $env:CUDAHOSTCXX
    $env:CXX = $env:CUDAHOSTCXX

    foreach ($tool in 'cl.exe', 'cmake.exe', 'ninja.exe', 'nvcc.exe', 'vcpkg.exe') {
        $resolved = Get-Command $tool -ErrorAction Stop
        Write-Host ('  {0,-12} {1}' -f $tool, $resolved.Source)
    }
    Write-Host "`nBuild environment ready (x64, Ninja, vcpkg, CUDA)." -ForegroundColor Green
    Write-Host 'Change to your project directory, then run:'
    Write-Host '  cmake -S . -B build'
    Write-Host '  cmake --build build'
}
