import numpy as np
import pandas as pd
from typing import Dict, Any

# Hardcoded data kept for fallback and validation
FALLBACK_TICKERS = [
    "UCB.BR", "EBAY", "BHP.AX", "7203.T", "PLTR", "SOF.BR", "005380.KS", "HOOD", "MSFT", "BMW.DE",
    "YZCFF", "NAB.AX", "BKNG", "ISP.MI", "MELI", "1810.HK", "AAPL", "VOW3.DE", "005930.KS", "6033.KL",
    "NVDA", "BBVA", "UNH", "BABA", "TSLA", "BBD", "SNAP", "SONY", "AMZN", "MC.PA",
    "KO", "HFG.DE", "TM", "PFE", "ORSTED.CO", "RACE", "INTC", "NTDOY", "AIR.PA", "JPM",
    "ENGI.PA", "000065.SZ", "NESN.SW", "FWONK", "PAG", "HMC", "MBG.DE", "NOK", "HPQ", "PINS",
    "9984.T", "ASML.AS", "MVST", "603080.SS", "ADBE", "NFLX", "SIE.DE", "SHEL", "CAT", "SAN.MC",
    "6752.T", "EA", "GE", "GLEN.L", "TCS.NS", "MQBKY", "BLK", "DIS", "OR.PA", "1211.HK",
    "TTE.PA", "F", "WMT", "META", "JNJ", "SHOP", "PEP", "BNP.PA", "000270.KS", "ADS.DE",
    "ZM", "LOTB.BR", "TNDM", "F9D.SI", "DIE.BR", "AMD", "NVO", "B4F.F", "KAEPY", "0992.HK",
    "QUBT", "BA.L", "ABBV", "ALIZF", "SAP", "V", "TATAMOTORS.NS", "ORCL", "MA", "CA.PA",
    "0700.HK", "2318.HK", "MAERSK-B.CO", "JPM", "MG", "COLR.BR", "7269.T", "SAF.PA", "JNJ", "IBE.MC",
    "FUJHY", "ANZ.AX", "NKE", "9601.T", "CS.PA", "UL", "HEIA.AS", "SBUX", "006400.KS", "SPOT",
    "0001.HK", "SMCI", "UPST", "CRSR", "BLNK", "PATH", "VALE", "EQNR.OL", "MMM", "BEP",
    "ADYEN.AS", "MT", "CEVA", "INFY", "PAH3.DE", "TDG", "DHR", "0857.HK", "ULVR.L", "2222.SR",
    "ROG.SW", "PG", "VOLV-B.ST", "6501.T", "TSM", "AI.PA", "PME.AX", "000333.SZ", "SOLB.BR", "BRK-B",
    "ACGBY", "BA", "IBM", "BP", "0941.HK", "E", "XTP.WA", "ELABS.OL", "KHC", "1821.T",
    "SAN.PA", "COST", "NHY.OL", "9983.T"
]



def load_market_data(file_path: str = "main techno.xlsx") -> Dict[str, Any]:
    print(f"Loading market data from {file_path}...")
    
    # Locate correct sheet (case insensitive)
    xl = pd.ExcelFile(file_path)
    sheet_name = None
    for name in xl.sheet_names:
        if "main" in name.lower().replace(" ", ""):
            sheet_name = name
            break
    if not sheet_name:
        raise ValueError("Could not find 'main' sheet in Excel file")
        
    df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)
    
    # Load Weights from Brouillon
    weights_sheet_name = None
    for name in xl.sheet_names:
        if "brouillon" in name.lower():
            weights_sheet_name = name
            break
    
    if weights_sheet_name:
        print(f"Loading weights0 from '{weights_sheet_name}'...")
        df_weights = pd.read_excel(file_path, sheet_name=weights_sheet_name, header=None)
        # We take fron sheet Brouillon, Column C (Index 2), Row 2 to 165 (Indices 1 to 164)
        weights0 = df_weights.iloc[1:165, 2].values.astype(float)
    else:
        print("WARNING: 'Brouillon' sheet not found. Using fallback/empty weights0.")
        weights0 = np.zeros(164)

    
    col_start = 2
    col_end = 166 # 164 items
    
    tickers = df.iloc[0, col_start:col_end].values.astype(str)
    mu = df.iloc[188, col_start:col_end].values.astype(float)
    beta = df.iloc[204, col_start:col_end].values.astype(float)
    alpha = df.iloc[212, col_start:col_end].values.astype(float)
    
    # Sanitize inputs (replace nans with 0 or handle them)
    weights0 = np.nan_to_num(weights0, nan=0.0)
    mu = np.nan_to_num(mu, nan=0.0)
    beta = np.nan_to_num(beta, nan=0.0)
    alpha = np.nan_to_num(alpha, nan=0.0)

    # Validate Ticker Order against Fallback
    # If loading UCB.BR, it should be first.
    if len(tickers) != len(FALLBACK_TICKERS):
        print(f"WARNING: Loacked tickers count ({len(tickers)}) != Fallback count ({len(FALLBACK_TICKERS)})")
    
    if tickers[0] != FALLBACK_TICKERS[0]:
        print(f"WARNING: First ticker mismatch! Loaded: {tickers[0]}, Expected: {FALLBACK_TICKERS[0]}")
    


    return {
        "tickers": np.array(tickers),
        "mu": mu,
        "weights0": weights0,
        "beta": beta,
        "alpha": alpha
    }
