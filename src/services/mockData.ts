// src/services/mockData.ts

export interface ChemicalData {
  id: number;
  name: string;
  smiles: string;
  score: number;  // Skor kecocokan/prediksi
  properties: string; // Sifat fisik (misal: Titik leleh, Kekuatan tarik)
  category: 'Thermoplastic' | 'Thermoset' | 'Elastomer' | 'Biodegradable';
  image?: string;
}

export const chemicalDatabase: ChemicalData[] = [
  {
    id: 1,
    name: "Polyethylene Terephthalate (PET)",
    smiles: "O=C(Oc1ccc(cc1)C(=O)OCCO)c2ccc(cc2)C(=O)O", // Contoh simplifikasi unit ulang
    score: 98.5,
    properties: "Melting Point: 260°C, High Tensile Strength",
    category: "Thermoplastic"
  },
  {
    id: 2,
    name: "Polylactic Acid (PLA)",
    smiles: "C(C(=O)O)C", // Unit dasar Lactic acid
    score: 95.2,
    properties: "Biodegradable, Melting Point: 150-160°C",
    category: "Biodegradable"
  },
  {
    id: 3,
    name: "Polypropylene (PP)",
    smiles: "CC(C)C", 
    score: 89.4,
    properties: "High Chemical Resistance, Toughness",
    category: "Thermoplastic"
  },
  {
    id: 4,
    name: "Polystyrene (PS)",
    smiles: "C(c1ccccc1)C",
    score: 82.1,
    properties: "Brittle, Transparent, Low Melting Point",
    category: "Thermoplastic"
  },
  {
    id: 5,
    name: "Polyurethane (PU)",
    smiles: "O=C(Nc1ccc(cc1)CC2=CC=C(C=C2)N=C=O)OCCO",
    score: 91.0,
    properties: "Flexible, Durable, Insulator",
    category: "Thermoset"
  }
];