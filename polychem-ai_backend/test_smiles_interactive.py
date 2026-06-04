#!/usr/bin/env python3
"""
Interactive SMILES Testing Script
Test backend predictions with various SMILES strings
"""

import requests
import json
import time
from typing import Tuple, List

BASE_URL = "http://127.0.0.1:8000"

# ANSI Colors
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

def print_header(text: str):
    print(f"\n{BLUE}{'='*80}{RESET}")
    print(f"{BLUE}{text:^80}{RESET}")
    print(f"{BLUE}{'='*80}{RESET}\n")

def print_success(text: str):
    print(f"{GREEN}✓ {text}{RESET}")

def print_error(text: str):
    print(f"{RED}✗ {text}{RESET}")

def test_smiles(smiles: str, description: str = "") -> bool:
    """Test a single SMILES string"""
    
    print(f"\n{BLUE}─────────────────────────────────────────────────────────────{RESET}")
    print(f"Testing: {smiles}")
    if description:
        print(f"Description: {description}")
    print(f"{BLUE}─────────────────────────────────────────────────────────────{RESET}")
    
    try:
        start = time.time()
        response = requests.post(
            f"{BASE_URL}/predict",
            json={"smiles": smiles},
            timeout=30,
            headers={"Content-Type": "application/json"}
        )
        elapsed = time.time() - start
        
        print(f"Response Time: {elapsed:.2f}s")
        print(f"Status Code: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            print_success("Prediction successful!")
            
            # New compound
            print(f"\n{YELLOW}New Compound Prediction:{RESET}")
            nc = data.get("new_compound", {})
            print(f"  Name: {nc.get('name', 'N/A')}")
            print(f"  Formula: {nc.get('formula', 'N/A')}")
            print(f"  Molecular Weight: {nc.get('molecular_weight', 0.0)} g/mol")
            print(f"  Tg (Glass Transition): {nc.get('tg', 0.0)}°C")
            print(f"  Justification: {nc.get('justifikasi', 'N/A')}")
            print(f"  Image: {nc.get('image_url', 'N/A')}")
            
            # Similar compounds
            print(f"\n{YELLOW}Similar Compounds (Top 3):{RESET}")
            similar = data.get("similar_compounds", [])
            for idx, comp in enumerate(similar[:3], 1):
                print(f"\n  #{idx} Similarity: {comp.get('similarity_percent', 0.0):.1f}%")
                print(f"     SMILES: {comp.get('smiles', 'N/A')}")
                print(f"     Name: {comp.get('name', 'N/A')}")
                print(f"     Tg: {comp.get('tg', 0.0)}°C")
                print(f"     Reason: {comp.get('justifikasi', 'N/A')}")
            
            return True
        else:
            error_msg = response.text[:200]
            print_error(f"API Error: {response.status_code}")
            print(f"Message: {error_msg}")
            return False
            
    except requests.exceptions.ConnectionError:
        print_error("Cannot connect to backend!")
        print(f"Make sure backend is running:")
        print(f"  python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000")
        return False
    except requests.exceptions.Timeout:
        print_error("Request timeout (>30 seconds)")
        return False
    except Exception as e:
        print_error(f"Error: {str(e)}")
        return False

def test_invalid_smiles(smiles: str, description: str = "") -> bool:
    """Test invalid SMILES (should return error)"""
    
    print(f"\n{BLUE}─────────────────────────────────────────────────────────────{RESET}")
    print(f"Testing Invalid SMILES: {smiles}")
    if description:
        print(f"Description: {description}")
    print(f"{BLUE}─────────────────────────────────────────────────────────────{RESET}")
    
    try:
        response = requests.post(
            f"{BASE_URL}/predict",
            json={"smiles": smiles},
            timeout=30
        )
        
        if response.status_code != 200:
            print_success("Correctly rejected invalid SMILES")
            print(f"Status: {response.status_code}")
            error_data = response.json()
            print(f"Error: {error_data.get('detail', 'N/A')}")
            return True
        else:
            print_error("Should have rejected invalid SMILES but didn't!")
            return False
            
    except Exception as e:
        print_error(f"Error: {str(e)}")
        return False

def main():
    """Main testing function"""
    
    print_header("POLYCHEM AI - SMILES PREDICTION TESTING")
    
    # Check if backend is running
    print("Checking backend health...\n")
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=5)
        if response.status_code == 200:
            print_success(f"Backend is running on {BASE_URL}")
            print_success(f"Health check: {response.json()}")
        else:
            print_error("Backend health check failed")
            return
    except Exception as e:
        print_error(f"Cannot reach backend: {str(e)}")
        return
    
    # Test valid SMILES
    print_header("VALID SMILES TESTS")
    
    valid_tests: List[Tuple[str, str]] = [
        ("C", "Methane - simplest organic"),
        ("CC", "Ethane - 2 carbons"),
        ("CCO", "Ethanol - common alcohol"),
        ("CCOC", "Ethyl methyl ether"),
        ("c1ccccc1", "Benzene - aromatic ring"),
        ("CC(=O)O", "Acetic acid - simple carboxylic acid"),
        ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin - well-known compound"),
        ("CCCCCCCCCCCCc1ccccc1", "Dodecylbenzene - long chain"),
    ]
    
    passed = 0
    for smiles, desc in valid_tests:
        if test_smiles(smiles, desc):
            passed += 1
    
    print(f"\n{BLUE}Valid SMILES Tests: {passed}/{len(valid_tests)} passed{RESET}")
    
    # Test invalid SMILES
    print_header("INVALID SMILES TESTS (Should be rejected)")
    
    invalid_tests: List[Tuple[str, str]] = [
        ("random string", "Not a valid SMILES"),
        ("###", "Invalid characters"),
        ("C@@@C", "Invalid bond notation"),
        ("", "Empty string"),
    ]
    
    invalid_passed = 0
    for smiles, desc in invalid_tests:
        # Skip empty string test (Pydantic min_length=1)
        if smiles and test_invalid_smiles(smiles, desc):
            invalid_passed += 1
    
    print(f"\n{BLUE}Invalid SMILES Tests: {invalid_passed}/{len(invalid_tests)-1} passed{RESET}")
    
    # Summary
    print_header("TEST SUMMARY")
    print(f"Valid SMILES: {passed}/{len(valid_tests)} passed")
    print(f"Invalid SMILES: {invalid_passed}/{len(invalid_tests)-1} passed")
    total = passed + invalid_passed
    total_tests = len(valid_tests) + len(invalid_tests) - 1
    print(f"\nTotal: {total}/{total_tests} tests passed")
    
    if total == total_tests:
        print_success("All tests passed!")
    else:
        print_error(f"{total_tests - total} tests failed")

if __name__ == "__main__":
    main()
