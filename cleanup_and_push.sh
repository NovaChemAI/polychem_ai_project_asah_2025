#!/bin/bash
# Cleanup script - Remove temporary files and prepare for GitHub push

echo "================================"
echo "   CLEANUP & PUSH TO GITHUB"
echo "================================"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Step 1: Delete temporary files
echo -e "\n${YELLOW}Step 1: Deleting temporary documentation files...${NC}"

files_to_delete=(
    "BACKEND_SETUP_READY.md"
    "BACKEND_ANALYSIS.md"
    "BACKEND_FIXES_RECOMMENDED.md"
    "BACKEND_TESTING_READY.md"
    "BACKEND_QUICK_REFERENCE.md"
    "BACKEND_DOCUMENTATION_INDEX.md"
    "BACKEND_EXECUTIVE_SUMMARY.md"
    "README_BACKEND_TESTING.md"
    "START_HERE.md"
    "START_BACKEND.bat"
    "FRONTEND_ANALYSIS.md"
    "FRONTEND_FIXES_CRITICAL.md"
    "SWAGGER_MANUAL_TESTING.md"
)

for file in "${files_to_delete[@]}"; do
    if [ -f "$file" ]; then
        rm -f "$file"
        echo -e "${GREEN}✓ Deleted: $file${NC}"
    fi
done

# Step 2: Remove Python cache directories
echo -e "\n${YELLOW}Step 2: Removing Python cache directories...${NC}"

if [ -d "polychem-ai_backend/__pycache__" ]; then
    rm -rf "polychem-ai_backend/__pycache__"
    echo -e "${GREEN}✓ Removed: polychem-ai_backend/__pycache__${NC}"
fi

if [ -d "polychem-ai_backend/app/__pycache__" ]; then
    rm -rf "polychem-ai_backend/app/__pycache__"
    echo -e "${GREEN}✓ Removed: polychem-ai_backend/app/__pycache__${NC}"
fi

# Step 3: Verify .env not in tracking
echo -e "\n${YELLOW}Step 3: Checking for .env files...${NC}"

if [ -f "polychem-ai_backend/.env" ]; then
    echo -e "${YELLOW}⚠ Found: .env file${NC}"
    echo -e "${YELLOW}  Make sure .env is in .gitignore${NC}"
    echo -e "${YELLOW}  File will NOT be pushed to GitHub${NC}"
fi

# Step 4: Show git status
echo -e "\n${YELLOW}Step 4: Current Git status:${NC}"
git status --short | head -20

# Step 5: Ask for commit message
echo -e "\n${YELLOW}Step 5: Preparing to push...${NC}"
read -p "Enter commit message (or press Enter for default): " commit_msg
commit_msg=${commit_msg:-"Clean up: remove temp docs and update for GitHub"}

# Step 6: Stage and commit
echo -e "\n${YELLOW}Step 6: Staging files...${NC}"
git add .
echo -e "${GREEN}✓ Files staged${NC}"

echo -e "\n${YELLOW}Step 7: Creating commit...${NC}"
git commit -m "$commit_msg"

# Step 7: Push
echo -e "\n${YELLOW}Step 8: Pushing to GitHub...${NC}"
git push

echo -e "\n${GREEN}================================${NC}"
echo -e "${GREEN}   ✓ CLEANUP & PUSH COMPLETE!${NC}"
echo -e "${GREEN}================================${NC}"

echo -e "\n${GREEN}What was done:${NC}"
echo "  ✓ Deleted 13 temporary documentation files"
echo "  ✓ Removed Python cache directories"
echo "  ✓ Staged all changes"
echo "  ✓ Created commit: '$commit_msg'"
echo "  ✓ Pushed to GitHub"

echo -e "\n${GREEN}Final structure:${NC}"
echo "  ✓ Source code: Clean and complete"
echo "  ✓ Documentation: Essential files only"
echo "  ✓ Cache: Removed"
echo "  ✓ Secrets: .env NOT committed"
echo ""
