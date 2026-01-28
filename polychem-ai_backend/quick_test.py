#!/usr/bin/env python3
"""
Quick SMILES Test - Minimal version
Just test one SMILES at a time
"""

import requests
import json

BASE_URL = "http://127.0.0.1:8000"

def test_smiles_simple(smiles: str):
    """Quick test a single SMILES"""
    
    print(f"\n{'='*60}")
    print(f"Testing SMILES: {smiles}")
    print(f"{'='*60}\n")
    
    try:
        # Send request
        response = requests.post(
            f"{BASE_URL}/predict",
            json={"smiles": smiles},
            timeout=30
        )
        
        print(f"Status: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            
            # Print new compound
            nc = data["new_compound"]
            print(f"\n✅ SUCCESS!\n")
            print(f"Name: {nc['name']}")
            print(f"Formula: {nc['formula']}")
            print(f"Tg: {nc['tg']}°C")
            print(f"Reason: {nc['justifikasi']}")
            
            # Print similar compounds
            print(f"\nTop 3 Similar:")
            for i, sc in enumerate(data["similar_compounds"][:3], 1):
                print(f"  {i}. {sc['name']} ({sc['similarity_percent']:.1f}% similar)")
                print(f"     {sc['justifikasi']}")
        else:
            print(f"\n❌ ERROR: {response.json().get('detail', 'Unknown error')}")
            
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    # Test some SMILES
    test_smiles_simple("C")      # Methane
    test_smiles_simple("CC")     # Ethane
    test_smiles_simple("CCO")    # Ethanol
    test_smiles_simple("c1ccccc1")  # Benzene
