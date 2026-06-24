#!/usr/bin/env python3
"""Quantitative value fidelity harness. Deterministic extraction, alignment, and
classification of quantitative values in patient-facing simplifications, with
severity-weighted metrics and Wilson confidence intervals."""
import json, re, math, csv

corpus = json.load(open("corpus.json"))["items"]
outputs = json.load(open("outputs.json"))["outputs"]

CATS = ["preserved","acceptably_simplified","distorted","unit_error","direction_error","omitted","fabricated"]

NUMRE = re.compile(r'(?<![\w.])(\d+(?:\.\d+)?)')

UNIT_PATTERNS = [
    ("mmol/L", r'mmol\s*/\s*l'),
    ("mg/dL",  r'mg\s*/\s*dl'),
    ("mcg",    r'microgram|mcg'),
    ("mmHg",   r'mmhg|mm\s*hg'),
    ("mg",     r'\bmg\b'),
    ("%",      r'%|percent'),
    ("C",      r'degree|celsius'),
    ("kg",     r'\bkg\b'),
    ("units",  r'\bunit'),
    ("weeks",  r'\bweek'),
    ("years",  r'\byear'),
    ("hours",  r'\bhour'),
    ("days",   r'\bday'),
]
DIR_WORDS = {"high":["high","higher","elevated","above","raised"],
             "low":["low","lower","below","reduced"]}

def extract(text):
    """Return list of dicts: {val, unit, ctx} for each number in text."""
    t = text.lower()
    out=[]
    for m in NUMRE.finditer(t):
        val=float(m.group(1))
        win=t[m.end():m.end()+14]
        unit=None
        for u,pat in UNIT_PATTERNS:
            if re.search(pat, win):
                # avoid tagging 'day' when it is 'times a day' frequency handled later
                unit=u; break
        out.append({"val":val,"unit":unit,"pos":m.start()})
    # shared-unit pass: in "X to Y unit" or "X and Y unit", a number lacking a
    # unit inherits the unit of a sibling number that appears just after it.
    for i,o in enumerate(out):
        if o["unit"] is None:
            seg=t[o["pos"]:o["pos"]+30]
            if re.search(r'\d+(?:\.\d+)?\s*(?:to|and|-|or)\s*\d', seg):
                for j in range(i+1, len(out)):
                    if out[j]["unit"] is not None and out[j]["pos"]-o["pos"]<30:
                        o["unit"]=out[j]["unit"]; break
    return out, t

def near(a,b,tol):
    return abs(a-b) <= tol + 1e-9

def freq_count(text):
    """Detect dosing frequency expressed in output (per day)."""
    t=text.lower()
    if "three times" in t or "3 times" in t: return 3
    if "twice" in t or "two times" in t or "2 times" in t: return 2
    if "once" in t or "one time" in t or "each morning" in t or "each night" in t or "a day" in t or "daily" in t or "each day" in t:
        return 1
    return None

def classify_item(item):
    out_text = outputs[item["id"]]["__cond__"]
    nums, t = extract(out_text)
    used=set()
    results=[]
    # percent reformulation: "X in 100" equals X percent
    in100 = re.search(r'(\d+(?:\.\d+)?)\s+in\s+100', t)
    in100_val = float(in100.group(1)) if in100 else None

    for g in item["gold"]:
        kind=g["kind"]; tol=g.get("tol",0); cat=None; note=""
        gn=g["nums"]
        if kind=="freq":
            fc=freq_count(out_text)
            if fc is None: cat="omitted"
            elif fc==gn[0]: cat="preserved"
            else: cat="distorted"; note=f"freq {fc} vs {gn[0]}"
            results.append((g,cat,note)); continue
        if kind in ("bp","ratio"):
            found=[any(near(o["val"],n,tol) for o in nums) for n in gn]
            if all(found): cat="preserved"
            elif any(found): cat="distorted"; note="partial multi-number"
            else: cat="omitted"
            results.append((g,cat,note)); continue

        n=gn[0]
        # candidate matches
        matched=[o for o in nums if near(o["val"],n,tol)]
        # percent reformulation handling
        if not matched and g["unit"]=="%" and in100_val is not None and near(in100_val,n,0):
            results.append((g,"acceptably_simplified","percent expressed as X in 100")); continue
        if matched:
            o=matched[0]
            # unit assessment
            gu=g["unit"]
            if gu and gu not in ("ratio",):
                if o["unit"]==gu:
                    cat="preserved"
                elif o["unit"] is None:
                    if g.get("unit_critical"): cat="unit_error"; note="critical unit dropped"
                    else: cat="acceptably_simplified"; note="unit dropped, value clear in context"
                else:
                    # different real unit
                    cat="unit_error"; note=f"unit {o['unit']} vs {gu}"
            else:
                cat="preserved"
            # direction check, scoped to a window around the matched value
            if g.get("direction"):
                opp = "low" if g["direction"]=="high" else "high"
                w0=max(0,o["pos"]-18); w1=o["pos"]+24
                window=t[w0:w1]
                if any(x in window for x in DIR_WORDS[opp]) and not any(x in window for x in DIR_WORDS[g["direction"]]):
                    cat="direction_error"; note="direction reversed"
            # rounding within tol but inexact -> acceptably simplified
            if cat=="preserved" and not any(abs(o["val"]-n)<1e-9 for o in matched):
                cat="acceptably_simplified"; note="rounded within tolerance"
            results.append((g,cat,note)); continue
        # no numeric match: distinguish distorted (concept present with wrong number) vs omitted
        concept=g["concept"].lower()
        concept_present = concept in t
        # a wrong number attached to concept => distorted
        if concept_present and nums:
            # heuristic: if a number sits within 25 chars of the concept and is not any gold value, call distorted
            ci=t.find(concept)
            close=[o for o in nums if abs(o["pos"]-ci)<25 and not any(near(o["val"],gn0,tol) for gn0 in gn)]
            if close:
                ratio=max(close[0]["val"]/n, n/close[0]["val"]) if close[0]["val"]>0 and n>0 else 0
                cat="distorted"; note=("order-of-magnitude " if ratio>=5 else "")+f"{close[0]['val']} vs {n}"
                results.append((g,cat,note)); continue
        cat="omitted"
        results.append((g,cat,note))

    # fabrication: output numbers not matching any gold and not part of percent reformulation
    gold_nums=[x for g in item["gold"] for x in g["nums"]]
    fabricated=[]
    for o in nums:
        if any(near(o["val"],gv,0.1) for gv in gold_nums): continue
        if in100_val is not None and o["val"]==100: continue   # 'in 100' scaffold
        # ignore years-of-age style? none here
        fabricated.append(o["val"])
    return results, fabricated

def wilson(k,n,z=1.96):
    if n==0: return (0,0,0)
    p=k/n
    den=1+z*z/n
    centre=(p+z*z/(2*n))/den
    half=(z*math.sqrt(p*(1-p)/n+z*z/(4*n*n)))/den
    return (p, max(0,centre-half), min(1,centre+half))

def run(cond):
    rows=[]; fab_total=0; out_nums_total=0
    for item in corpus:
        outputs[item["id"]]["__cond__"]=outputs[item["id"]][cond]
        res, fab = classify_item(item)
        nums,_=extract(outputs[item["id"]][cond])
        out_nums_total+=len(nums); fab_total+=len(fab)
        for g,cat,note in res:
            rows.append({"item":item["id"],"family":item["family"],"gid":g["gid"],
                         "concept":g["concept"],"gold":"/".join(map(str,g["nums"]))+(" "+(g["unit"] or "") if g["unit"] else ""),
                         "severity":g["severity"],"category":cat,"note":note,
                         "fabricated_in_item":";".join(map(str,fab)) if fab else ""})
    return rows, fab_total, out_nums_total

def summarise(rows, fab_total, out_nums_total):
    N=len(rows)
    sev_total=sum(r["severity"] for r in rows)
    def cnt(c): return sum(1 for r in rows if r["category"]==c)
    preserved=cnt("preserved")+cnt("acceptably_simplified")
    distorted=cnt("distorted")+cnt("unit_error")+cnt("direction_error")
    omitted=cnt("omitted")
    # critical subset (severity 3)
    crit=[r for r in rows if r["severity"]==3]
    crit_err=sum(1 for r in crit if r["category"] in ("distorted","unit_error","direction_error","omitted"))
    # severity-weighted fidelity index
    good_w=sum(r["severity"] for r in rows if r["category"] in ("preserved","acceptably_simplified"))
    fab_penalty=fab_total*2  # each fabrication penalised at moderate weight
    qfi=max(0,(good_w - fab_penalty))/sev_total if sev_total else 0
    return {
      "n_values":N, "preserved":preserved, "distorted":distorted, "omitted":omitted,
      "fabricated":fab_total, "out_numbers":out_nums_total, "n_critical":len(crit),
      "crit_err":crit_err,
      "preservation_rate":wilson(preserved,N),
      "distortion_rate":wilson(distorted,N),
      "omission_rate":wilson(omitted,N),
      "fabrication_rate":wilson(fab_total,out_nums_total),
      "critical_value_error_rate":wilson(crit_err,len(crit)),
      "qfi":qfi,
      "by_cat":{c:cnt(c) for c in CATS}
    }

allrows={}
summary={}
for cond in ["plain","preserve"]:
    rows, fab, outn = run(cond)
    allrows[cond]=rows
    summary[cond]=summarise(rows, fab, outn)
    with open(f"audit_{cond}.csv","w",newline="") as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

json.dump(summary, open("results.json","w"), indent=2)

# pretty print
for cond in ["plain","preserve"]:
    s=summary[cond]
    print(f"\n===== CONDITION: {cond} =====")
    print(f"values scored: {s['n_values']}  | output numbers: {s['out_numbers']}  | critical values: {s['n_critical']}")
    print("category counts:", s["by_cat"])
    for k in ["preservation_rate","distortion_rate","omission_rate","fabrication_rate","critical_value_error_rate"]:
        p,lo,hi=s[k]; print(f"  {k:28s}: {p*100:5.1f}%  (95% CI {lo*100:.1f} to {hi*100:.1f})")
    print(f"  Quantitative Fidelity Index : {s['qfi']:.3f}")
