"""R2 [ПОМ] — результаты проверки тетрад, сохранённые для [ОСН]."""
import json

res = {
    "TITLE": "R2 [ПОМ]: torsion scalar хорошей vs диагональной тетрады (Morris-Thorne, A=1, B=1/sqrt(1-b0^2/r^2))",
    "closed_form_good_tetrad": {
        "T_good(r)": "4*(r - sqrt(r^2 - b0^2))^2 / r^4",
        "asymptotics": "~ b0^4/r^6  (r>>b0)",
        "verified_numeric": True,
        "angular_independence": "Confirmed numerically at several th,ph (scalar invariant)",
    },
    "diagonal_tetrad_R1": {
        "T_diag(r)": "1/r^2",
        "b0_dependence": "NONE — T_diag does NOT depend on b0/shape function",
        "conclusion": "BAD TETRAD: torsion scalar does not encode b(r); f(T)-dynamics built on diagonal tetrad cannot test the wormhole throat. Prescribe recomputing R1 with a good (symmetric) tetrad.",
        "literature": "Tamanini & Boehmer, arXiv:1204.4593; Boehmer-Harko-Lobo PRD85 044033 the good-tetrad approach",
    },
    "numeric_refs_good_tetrad_b0_0_1": {
        "r=0.2": 1.794919, "r=0.3": 0.1453691, "r=0.5": 0.006531292,
        "r=1.0": 0.0001005031, "r=2.0": 1.564456e-06, "r=3.0": 1.372505e-07,
    },
    "next": "R2 [ПОМ]: NEC/WEC scan for f(T)=T+alpha*T^2 and power-law f(T)=T+lambda*(-T)^n on good tetrad",
}
with open("ft2_tetrad_findings.json", "w") as f:
    json.dump(res, f, indent=2)
print("saved ft2_tetrad_findings.json")
print("T_good(r) = 4*(r - sqrt(r^2 - b0^2))^2 / r^4")
print("T_diag(r) = 1/r^2  (b0-independent => bad tetrad, throat invisible to f(T))")