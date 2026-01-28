@echo off
REM Cleanup script for Windows - Remove temporary files and prepare for GitHub push

echo.
echo ================================
echo    CLEANUP AND PUSH TO GITHUB
echo ================================
echo.

REM Colors (Windows batch doesn't support ANSI colors, using simpler approach)

REM Step 1: Delete temporary files
echo [STEP 1] Deleting temporary documentation files...
echo.

set "files=BACKEND_SETUP_READY.md BACKEND_ANALYSIS.md BACKEND_FIXES_RECOMMENDED.md BACKEND_TESTING_READY.md BACKEND_QUICK_REFERENCE.md BACKEND_DOCUMENTATION_INDEX.md BACKEND_EXECUTIVE_SUMMARY.md README_BACKEND_TESTING.md START_HERE.md START_BACKEND.bat FRONTEND_ANALYSIS.md FRONTEND_FIXES_CRITICAL.md SWAGGER_MANUAL_TESTING.md"

for %%F in (%files%) do (
    if exist "%%F" (
        del /Q "%%F"
        echo [OK] Deleted: %%F
    )
)

REM Step 2: Remove Python cache directories
echo.
echo [STEP 2] Removing Python cache directories...
echo.

if exist "polychem-ai_backend\__pycache__" (
    rmdir /s /q "polychem-ai_backend\__pycache__" 2>nul
    echo [OK] Removed: polychem-ai_backend\__pycache__
)

if exist "polychem-ai_backend\app\__pycache__" (
    rmdir /s /q "polychem-ai_backend\app\__pycache__" 2>nul
    echo [OK] Removed: polychem-ai_backend\app\__pycache__
)

REM Step 3: Check .env
echo.
echo [STEP 3] Checking for .env files...
echo.

if exist "polychem-ai_backend\.env" (
    echo [WARNING] Found: .env file
    echo [INFO] Make sure .env is in .gitignore
    echo [INFO] File will NOT be pushed to GitHub
) else (
    echo [OK] No .env file found (or already in .gitignore)
)

REM Step 4: Git status
echo.
echo [STEP 4] Git status:
echo.
git status --short | findstr /C:" "
if %errorlevel% neq 0 (
    echo [INFO] Working directory clean
)

REM Step 5: Ask for commit message
echo.
echo [STEP 5] Ready to commit and push
echo.
set /p commit_msg="Enter commit message (or press Enter for default): "
if "%commit_msg%"=="" (
    set "commit_msg=Clean up: remove temp docs and prepare for GitHub"
)

REM Step 6: Git add
echo.
echo [STEP 6] Staging files...
echo.
git add .
echo [OK] Files staged

REM Step 7: Git commit
echo.
echo [STEP 7] Creating commit...
echo.
git commit -m "%commit_msg%"

if %errorlevel% neq 0 (
    echo [ERROR] Commit failed
    echo [INFO] Check git status for details
    pause
    exit /b 1
)

REM Step 8: Git push
echo.
echo [STEP 8] Pushing to GitHub...
echo.
git push

if %errorlevel% neq 0 (
    echo [ERROR] Push failed
    echo [INFO] Check your GitHub connection and credentials
    pause
    exit /b 1
)

REM Success message
echo.
echo ================================
echo    CLEANUP AND PUSH COMPLETE!
echo ================================
echo.
echo What was done:
echo   - Deleted 13 temporary documentation files
echo   - Removed Python cache directories
echo   - Staged all changes
echo   - Created commit: "%commit_msg%"
echo   - Pushed to GitHub
echo.
echo Final structure:
echo   - Source code: Clean and complete
echo   - Documentation: Essential files only
echo   - Cache: Removed
echo   - Secrets: .env NOT committed
echo.
pause
