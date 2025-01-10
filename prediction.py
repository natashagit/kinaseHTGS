import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
import numpy as np
from warnings import filterwarnings

filterwarnings('ignore')

def check_solubility(dataset, props):

        dataset = "filtered_SOL.csv"  
        data = pd.read_csv(dataset)

        X = data.iloc[:, :-1] 
        y = data.iloc[:, -1]   # label


        weights = [0.25, 0.125, 0.125, 0.5] 
        X_weighted = X * weights  

        # for regressor
        X_train, X_test, y_train, y_test = train_test_split(X_weighted, y, test_size=0.2, random_state=42)

        model = RandomForestRegressor(random_state=42)
        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)

        r2 = r2_score(y_test, y_pred)
        print(f"R-squared on Test Set: {r2:.2f}")

        props = np.array(props).reshape(1, -1)
        weighted_props = props * weights
        sol = model.predict(weighted_props)
        print(f"solubility: {sol}")

        if sol > 0:
            label = "Highly Soluble"
        elif -2.303 < sol <= 0:
            label = "Soluble"
        elif -4.605 < sol <= -2.303:
            label =  "Low Solubility"
        else:
            label =  "Insoluble"
        
        return label



