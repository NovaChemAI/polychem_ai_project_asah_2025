from pydantic import BaseModel, Field
from typing import List


class RecommendRequest(BaseModel):
    smiles: str = Field(..., min_length=1)


class NewCompoundOut(BaseModel):
    name: str
    smiles: str
    formula: str = ""
    molecular_weight: float = 0.0
    tg_justification: str = ""  
    tg: float = 0.0
    pid: str = ""
    polymer_class: str = ""
    justifikasi: str
    fingerprint_length: int
    image_filename: str
    image_url: str


class SimilarCompoundOut(BaseModel):
    rank: int
    smiles: str
    name: str = ""
    formula: str = ""
    molecular_weight: float = 0.0
    tg: float = 0.0
    pid: str = ""
    polymer_class: str = ""
    similarity_score: float
    similarity_percent: float = 0.0
    justifikasi: str
    image_filename: str
    image_url: str


class RecommendResponse(BaseModel):
    status: str
    input_smiles: str
    new_compound: NewCompoundOut
    similar_compounds: List[SimilarCompoundOut]


class HistoryItemOut(BaseModel):
    input_smiles: str
    result: RecommendResponse