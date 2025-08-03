@echo off

:: Check if at least one argument (the source file) is provided
if "%~1"=="" (
    echo Usage: %~nx0 ^<source_file.f90^> [output_executable_name]
    exit /b 1
)

:: Use the provided source file and output executable name or default to "a.exe"
set "SOURCE_FILE=%~1"
set "OUTPUT_NAME=%~2"

if "%OUTPUT_NAME%"=="" (
    set "OUTPUT_NAME=a.exe"
)

:: Run the compilation command
ifx -o "%OUTPUT_NAME%" "%SOURCE_FILE%" -I../LibDualzn128 ../LibDualzn128/libdualzn.lib

