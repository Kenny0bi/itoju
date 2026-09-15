"""Plain names for traits and FinnGen endpoints, and the order figures show them in.

One place, so every figure and table names a definition the same way.
"""

TRAIT_LABEL = {
    "ASD": "Autism",
    "SCZ": "Schizophrenia",
    "BIP": "Bipolar",
    "MDD": "Depression",
    "PTSD": "PTSD",
    "MDD_EHR": "Depression, EHR",
    "MDD_Clin": "Depression, clinical",
}
TRAIT_SHORT = {"ASD": "ASD", "SCZ": "SCZ", "BIP": "BIP", "MDD": "MDD", "PTSD": "PTSD",
               "MDD_EHR": "MDD-EHR", "MDD_Clin": "MDD-Clin"}
TRAIT_ORDER = ["ASD", "SCZ", "BIP", "MDD", "PTSD", "MDD_EHR", "MDD_Clin"]

# condition -> [(endpoint, label)], definition sets first, neighbouring conditions after
GROUPS = [
    ("Epilepsy", [
        ("G6_EPLEPSY", "Any epilepsy"),
        ("FE", "Focal"),
        ("FE_STRICT", "Focal, strict"),
        ("FE_MODE", "Focal, mode"),
        ("GE", "Generalized"),
        ("GE_STRICT", "Generalized, strict"),
        ("GE_MODE", "Generalized, mode"),
    ]),
    ("Sleep apnoea", [
        ("G6_SLEEPAPNO", "Hospital records"),
        ("G6_SLEEPAPNO_INCLAVO", "+ primary care"),
        ("SLEEP", "Any sleep disorder"),
    ]),
    ("Insomnia", [
        ("F5_INSOMNIA", "Insomnia (F51.0, G47.0)"),
        ("KRA_PSY_SLEEP_NONORG_EXMORE", "Nonorganic sleep (F51)"),
    ]),
    ("Constipation", [
        ("K11_CONSTIPATION", "Constipation or laxatives"),
        ("K11_OTHFUNC", "Functional bowel (K59)"),
    ]),
    ("ADHD", [
        ("F5_ADHD", "F90.0"),
        ("KRA_PSY_HYPERKIN_EXMORE", "F90, strict controls"),
    ]),
    ("Intellectual disability", [
        ("F5_MILDRET", "Mild (F70)"),
        ("KRA_PSY_MENTALRET_EXMORE", "Any (F7), strict controls"),
    ]),
    ("Neighbouring conditions", [
        ("G6_STATUSEPI", "Status epilepticus"),
        ("G6_SLEEPDISOTH", "Other sleep disorders"),
        ("F5_SLEEP_NOS", "Sleep disorder, unspecified"),
        ("K11_IBS", "Irritable bowel"),
        ("K11_FUNCDYSP", "Functional dyspepsia"),
        ("K11_REFLUX", "Reflux"),
        ("N14_NEUROMUSCDYSBLADD", "Neurogenic bladder (N31)"),
        ("N14_OTHBLADD", "Other bladder (N32)"),
        ("KRA_PSY_DEVWIDE_EXMORE", "Autism spectrum (F84)"),
        ("KRA_PSY_AUTISM_EXMORE", "Autism (F84.0, F84.5)"),
    ]),
]
ENDPOINT_LABEL = {e: lab for _, items in GROUPS for e, lab in items}
ENDPOINT_GROUP = {e: g for g, items in GROUPS for e, _ in items}
ENDPOINT_ORDER = [e for _, items in GROUPS for e, _ in items]
