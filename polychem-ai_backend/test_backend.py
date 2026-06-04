#!/usr/bin/env python3
"""
Test Script untuk Backend PolyChem AI
Testing semua endpoints dan functionality
"""

import requests
import json
import time
from typing import Dict, Any

BASE_URL = "http://127.0.0.1:8000"

# ANSI Colors untuk output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

def print_result(test_name: str, success: bool, message: str = ""):
    status = f"{GREEN}✓ PASS{RESET}" if success else f"{RED}✗ FAIL{RESET}"
    print(f"{status} | {test_name}")
    if message:
        print(f"       └─ {message}")

def print_section(title: str):
    print(f"\n{BLUE}{'='*60}{RESET}")
    print(f"{BLUE}{title:^60}{RESET}")
    print(f"{BLUE}{'='*60}{RESET}\n")

# ============================================================
# TEST 1: Health Check
# ============================================================
def test_health():
    """Test endpoint /health"""
    print_section("TEST 1: Health Check")
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=5)
        success = response.status_code == 200
        data = response.json()
        print_result("Health Check", success, f"Status: {response.status_code}, Response: {data}")
        return success
    except Exception as e:
        print_result("Health Check", False, str(e))
        return False

# ============================================================
# TEST 2: Predict dengan SMILES Valid
# ============================================================
def test_predict_valid_smiles():
    """Test /predict dengan SMILES yang valid"""
    print_section("TEST 2: Predict dengan Valid SMILES")
    
    # Test dengan beberapa SMILES sederhana
    test_smiles = [
        ("C", "Methane - smallest organic"),
        ("CC", "Ethane"),
        ("CCO", "Ethanol"),
        ("c1ccccc1", "Benzene"),
    ]
    
    all_passed = True
    for smiles, description in test_smiles:
        try:
            print(f"\n  Testing: {description} ({smiles})")
            payload = {"smiles": smiles}
            response = requests.post(
                f"{BASE_URL}/predict",
                json=payload,
                timeout=30
            )
            
            success = response.status_code == 200
            if success:
                data = response.json()
                print_result(f"Predict {description}", True, 
                           f"Status: {response.status_code}")
                # Print hasil prediksi
                if "new_compound" in data:
                    nc = data["new_compound"]
                    print(f"       ├─ Name: {nc.get('name', 'N/A')}")
                    print(f"       ├─ Formula: {nc.get('formula', 'N/A')}")
                    print(f"       ├─ Tg: {nc.get('tg', 'N/A')}")
                    print(f"       └─ Similar compounds: {len(data.get('similar_compounds', []))}")
            else:
                error_msg = response.text[:100]
                print_result(f"Predict {description}", False, 
                           f"Status: {response.status_code}, Error: {error_msg}")
                all_passed = False
                
        except requests.exceptions.Timeout:
            print_result(f"Predict {description}", False, "Request timeout (>30s)")
            all_passed = False
        except Exception as e:
            print_result(f"Predict {description}", False, str(e))
            all_passed = False
    
    return all_passed

# ============================================================
# TEST 3: Predict dengan SMILES Invalid
# ============================================================
def test_predict_invalid_smiles():
    """Test /predict dengan SMILES yang invalid"""
    print_section("TEST 3: Predict dengan Invalid SMILES")
    
    invalid_smiles = [
        ("INVALID123", "Random string"),
        ("", "Empty string"),
        ("XxXxXx", "Invalid characters"),
    ]
    
    all_passed = True
    for smiles, description in invalid_smiles:
        try:
            print(f"\n  Testing: {description} ({smiles})")
            payload = {"smiles": smiles}
            response = requests.post(
                f"{BASE_URL}/predict",
                json=payload,
                timeout=10
            )
            
            # Invalid SMILES seharusnya return error (4xx atau 5xx)
            success = response.status_code >= 400
            print_result(f"Handle Invalid: {description}", success,
                       f"Status: {response.status_code}")
            if not success:
                all_passed = False
                
        except Exception as e:
            # Exception juga diterima untuk invalid input
            print_result(f"Handle Invalid: {description}", True, f"Exception caught: {str(e)[:50]}")
    
    return all_passed

# ============================================================
# TEST 4: History Endpoint
# ============================================================
def test_history():
    """Test endpoint /history"""
    print_section("TEST 4: History Endpoint")
    try:
        response = requests.get(f"{BASE_URL}/history", timeout=5)
        success = response.status_code == 200
        data = response.json()
        
        is_list = isinstance(data, list)
        print_result("Get History", success and is_list, 
                   f"Status: {response.status_code}, Items: {len(data) if is_list else 'N/A'}")
        
        if is_list and len(data) > 0:
            print(f"       └─ Recent item keys: {list(data[0].keys())}")
        
        return success and is_list
    except Exception as e:
        print_result("Get History", False, str(e))
        return False

# ============================================================
# TEST 5: Response Format Validation
# ============================================================
def test_response_format():
    """Test bahwa response format sesuai dengan schema"""
    print_section("TEST 5: Response Format Validation")
    
    try:
        payload = {"smiles": "CCO"}
        response = requests.post(f"{BASE_URL}/predict", json=payload, timeout=30)
        
        if response.status_code != 200:
            print_result("Response Format", False, f"Non-200 status: {response.status_code}")
            return False
        
        data = response.json()
        
        # Validasi struktur response
        required_fields = ["status", "input_smiles", "new_compound", "similar_compounds"]
        has_all_fields = all(field in data for field in required_fields)
        
        print_result("Has required fields", has_all_fields, 
                   f"Fields: {required_fields}")
        
        if has_all_fields:
            # Validasi new_compound structure
            nc = data["new_compound"]
            nc_fields = ["name", "smiles", "formula", "molecular_weight", "tg", "justifikasi", "image_url"]
            has_nc_fields = all(field in nc for field in nc_fields)
            
            print_result("NewCompound has required fields", has_nc_fields,
                       f"Fields: {nc_fields}")
            
            # Validasi similar_compounds structure
            sim = data["similar_compounds"]
            if isinstance(sim, list) and len(sim) > 0:
                sim_fields = ["rank", "smiles", "name", "similarity_score", "similarity_percent", "justifikasi", "image_url"]
                has_sim_fields = all(field in sim[0] for field in sim_fields)
                
                print_result("SimilarCompound has required fields", has_sim_fields,
                           f"Fields: {sim_fields}")
                
                return has_all_fields and has_nc_fields and has_sim_fields
        
        return has_all_fields
        
    except Exception as e:
        print_result("Response Format", False, str(e))
        return False

# ============================================================
# TEST 6: Performance Test
# ============================================================
def test_performance():
    """Test kecepatan response"""
    print_section("TEST 6: Performance Test")
    
    try:
        payload = {"smiles": "CCO"}
        
        start = time.time()
        response = requests.post(f"{BASE_URL}/predict", json=payload, timeout=60)
        elapsed = time.time() - start
        
        success = response.status_code == 200
        is_fast = elapsed < 30  # Harus < 30 detik
        
        print_result(f"Predict Performance", success and is_fast,
                   f"Response time: {elapsed:.2f}s (target: <30s)")
        
        return success and is_fast
        
    except Exception as e:
        print_result("Performance Test", False, str(e))
        return False

# ============================================================
# TEST 7: Concurrent Requests (Simple)
# ============================================================
def test_concurrent():
    """Test multiple requests"""
    print_section("TEST 7: Multiple Sequential Requests")
    
    smiles_list = ["C", "CC", "CCO"]
    all_passed = True
    
    for i, smiles in enumerate(smiles_list, 1):
        try:
            response = requests.post(
                f"{BASE_URL}/predict",
                json={"smiles": smiles},
                timeout=30
            )
            success = response.status_code == 200
            print_result(f"Request {i} ({smiles})", success,
                       f"Status: {response.status_code}")
            if not success:
                all_passed = False
        except Exception as e:
            print_result(f"Request {i} ({smiles})", False, str(e))
            all_passed = False
    
    return all_passed

# ============================================================
# TEST 8: Edge Cases
# ============================================================
def test_edge_cases():
    """Test edge cases"""
    print_section("TEST 8: Edge Cases")
    
    all_passed = True
    
    # Test dengan SMILES dengan spasi
    try:
        response = requests.post(
            f"{BASE_URL}/predict",
            json={"smiles": "  CCO  "},
            timeout=30
        )
        success = response.status_code == 200
        print_result("SMILES with whitespace", success)
        if not success:
            all_passed = False
    except Exception as e:
        print_result("SMILES with whitespace", False, str(e))
        all_passed = False
    
    # Test dengan SMILES case sensitivity
    try:
        response = requests.post(
            f"{BASE_URL}/predict",
            json={"smiles": "c1ccccc1"},
            timeout=30
        )
        success = response.status_code == 200
        print_result("Aromatic SMILES (lowercase c)", success)
        if not success:
            all_passed = False
    except Exception as e:
        print_result("Aromatic SMILES", False, str(e))
        all_passed = False
    
    return all_passed

# ============================================================
# Main Execution
# ============================================================
def main():
    print(f"\n{BLUE}{'*'*60}{RESET}")
    print(f"{BLUE}{'PolyChem AI Backend Test Suite':^60}{RESET}")
    print(f"{BLUE}{'*'*60}{RESET}")
    print(f"Target URL: {BASE_URL}\n")
    
    print(f"{YELLOW}Waiting 2 seconds for backend to stabilize...{RESET}")
    time.sleep(2)
    
    results = {}
    
    # Run all tests
    results["Health Check"] = test_health()
    results["Valid SMILES"] = test_predict_valid_smiles()
    results["Invalid SMILES"] = test_predict_invalid_smiles()
    results["History"] = test_history()
    results["Response Format"] = test_response_format()
    results["Performance"] = test_performance()
    results["Concurrent"] = test_concurrent()
    results["Edge Cases"] = test_edge_cases()
    
    # Summary
    print_section("TEST SUMMARY")
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for test_name, result in results.items():
        status = f"{GREEN}PASS{RESET}" if result else f"{RED}FAIL{RESET}"
        print(f"{status} | {test_name}")
    
    print(f"\n{BLUE}{'='*60}{RESET}")
    print(f"Total: {passed}/{total} tests passed")
    
    if passed == total:
        print(f"{GREEN}All tests passed!{RESET}")
    else:
        print(f"{RED}{total - passed} test(s) failed{RESET}")
    
    print(f"{BLUE}{'='*60}{RESET}\n")

if __name__ == "__main__":
    main()
