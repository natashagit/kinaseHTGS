from flask import Flask, render_template, request
from chemistry import compute_all

app = Flask(__name__)


# Take the input string from user
# Get the Property values: Polar Surface Area(PSA), Number of Rotations(Nrot), Hydrogen Bond Donor(HBD), Hydrogen Bond Acceptor(HBA) 
@app.route("/", methods=["GET"])
def get_properties():
    smiles_query = request.args.get("query", "")
    if smiles_query!="":
        properties = compute_all(smiles_query)
    return render_template(
        "base.html", smiles_query=smiles_query, Nrot_value = properties['nrot'], HBD_value = properties['hbd'], HBA_value = properties['hba']
    )


# Calculate Solubility of molecule in water
# Render the property values, solubility and image for visualization into html



# @app.route("/", methods=["GET"])
# def root():
#     return render_template("base.html")


# Thresholds

def solubility(props):

    mw, psa, n_rot, hbd, hba = props['MW'], props['psa'], props['nrot'], props['hbd'], props['hba']
    PSA, NROT, HBD, HBA, A_D_ratio = False, False, False, False, False

    if mw <= 400:
    
        if psa > 0.2*mw:
            print ('psa within range')
            PSA = True
        else:
            print ('low psa')
        
        if hba > hbd:
            print ('HBA > HBD')
            A_D_ratio = True
        else:
            print ('HBD higher than HBA')

        if hbA < 10 and hbA > 7:
            print('hbA within range')
            HBA = True
        elif hbA < 7:
            print ('hba low. probably not')
        else:
            print ('hba > 10. could be yes')
            HBA = True
            
        if hbd < 3:
            print ('hbD within range')
            HBD = True
        
        if n_rot<6:
            print ('n_rot is within the threshold')
            NROT = True

    
    elif MW > 400 and MW < 550:

        if psa > 0.2*mw:
            print ('psa within range')
            PSA = True
        else:
            print ('low psa')
        
        if hba > hbd:
            print ('HBA > HBD')
            A_D_ratio = True
        else:
            print ('HBD higher than HBA')

        if hbA < 8 and hbA > 6:
            print('hbA within range')
            HBA = True
        elif hbA < 6:
            print ('hba low. probably not')
        else:
            print ('hba > 8. could be yes')
            HBA = True
            
        if hbd < 3:
            print ('hbD within range')
            HBD = True
        
        if n_rot<9:
            print ('n_rot is within the threshold')
            NROT = True

    else:
            if psa > 0.2*mw:
            print ('psa within range')
            PSA = True
        else:
            print ('low psa')
        
        if hba > hbd:
            print ('HBA > HBD')
            A_D_ratio = True
        else:
            print ('HBD higher than HBA')

        if hbA < 8 and hbA > 6:
            print('hbA within range')
            HBA = True
        elif hbA < 6:
            print ('hba low. probably not')
        else:
            print ('hba > 8. could be yes')
            HBA = True
            
        if hbd < 3:
            print ('hbD within range')
            HBD = True
        
        if n_rot<9:
            print ('n_rot is within the threshold')
            NROT = True

        return True