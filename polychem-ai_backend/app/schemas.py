from pydantic import BaseModel, Field
from typing import List, Dict, Any


class RecommendRequest(BaseModel):
    smiles: str = Field(..., min_length=1)


class NewCompoundOut(BaseModel):
    name: str
    smiles: str
    formula: str
    molecular_weight: float
    tg: float
    pid: str
    polymer_class: str
    justifikasi: str
    fingerprint_length: int
    image_filename: str
    image_url: str


class SimilarCompoundOut(BaseModel):
    rank: int
    smiles: str
    name: str
    formula: str
    molecular_weight: float
    tg: float
    pid: str
    polymer_class: str
    similarity_score: float
    similarity_percent: float
    justifikasi: str
    image_filename: str
    image_url: str


class RecommendResponse(BaseModel):
    status: str
    input_smiles: str
    new_compound: NewCompoundOut
    similar_compounds: List[SimilarCompoundOut]


# ✅ untuk endpoint /history
class HistoryItemOut(BaseModel):
    input_smiles: str
    result: RecommendResponse
