# Desenvolvido por:
#   - Leonardo Rossi Leão
#   - Matheus Meirelles Onofre Martins

from matplotlib import pyplot as plt
from scipy import optimize as op
import numpy as np

# Audi A4 1.8

def simulacao(tamb, pamb, op="", printTable=0):

    if op == "potencia":
        params = [0.866288841, 7.06979930, 952.584869, 1.15104888, 1.36923106, 1.25428223, 0.920010780, 0.899922327, 2458.76386]
    else:
        params = [0.9, 8.54308881, 939.255822, 1.15364973, 1.32752461, 1.24861101, 0.97, 0.803650162, 2414.5877]

    pdc, deltaT, tres, presp1, compexo, expexpo, fii, fii_gas, digitar_c = params

    referencia = {
        "potencia_efetiva": 56.1/1.014,        # HP
        "torque": 80.41453                   # N*m
    }

    PCI = {
        "gasolina": 10.4*4.1868,              # MJ/kg
        "etanol": 6.75*4.1868                 # MJ/kg
    }

    num_cilindros = 4
    diametro = 0.076                   # m
    curso = 0.0548                      # m
    tempo = 4
    x = tempo/2
    taxa_compressao = 9.5
    if op == "potencia":
        rotacao = 6000
    else:
        rotacao = 3250
    cilindrada = 994*1e-6

    fuel = {"tipo": "etanol", "alfa": 0.9}
    h2_co = 0.45

    if fuel["tipo"] == "gasolina":
        fuel["substancia"] = {"C": 8, "H": 18, "O": 0}
        fuel["peso"] = 114
    elif fuel["tipo"] == "etanol":
        fuel["substancia"] = {"C": 2, "H": 6, "O": 1}
        fuel["peso"] = 46

    # %%        Passo 1

    # Combustão estequiométrica
    comb_esteq = {
        "CO2": fuel["substancia"]["C"],
        "H2O": fuel["substancia"]["H"]/2
    }
    comb_esteq["O2"] = (2*comb_esteq["CO2"] + comb_esteq["H2O"] - fuel["substancia"]["O"])/2

    # Coeficientes da combustão real
    comb_coef = {
        "fuel": 1,
        "O2": comb_esteq["O2"]*fuel["alfa"],
        "N2": comb_esteq["O2"]*fuel["alfa"]*3.76,
        "H2": (2*fuel["alfa"]*comb_esteq["O2"] - fuel["substancia"]["H"]/2 - 2*fuel["substancia"]["C"])/(-2 + 1 - 1/h2_co)
    }
    comb_coef["H2O"] =  fuel["substancia"]["H"]/2 - comb_coef["H2"]
    comb_coef["CO"] = comb_coef["H2"]/h2_co
    comb_coef["CO2"] = fuel["substancia"]["C"] - comb_coef["CO"]

    comb_coef["produtos"] = (comb_coef["N2"] + comb_coef["CO2"] + comb_coef["CO"] + comb_coef["H2O"] + comb_coef["H2"])/fuel["peso"]    
    comb_coef["reagentes"] = (comb_coef["fuel"] + comb_coef["O2"] + comb_coef["N2"])/fuel["peso"]
    comb_coef["ar"] = (comb_coef["O2"] + comb_coef["N2"])/fuel["peso"]


    # %%        Passo 1: Admissao

    admissao = {
        "Tamb": tamb,           # Temperatura ambiente
        "pamb": pamb,         # Pressão ambiente
        "pdc": pdc,             # Perda de carga, varia entre: 0.8 a 0.9
        "deltaT": deltaT,            # Variacao de temperatura, entre 0 a 20ºC
        "tres": tres,            # , entre 900 e 1000 K
        "pres/p1": presp1         # , entre 1.1 a 1.25
    }

    resultados = {"p1": admissao["pamb"]*admissao["pdc"]}       
    admissao["pres"] = resultados["p1"]*admissao["pres/p1"]
    admissao["epson"] = taxa_compressao

    resultados["efe"] = ((admissao["Tamb"] + admissao["deltaT"])/admissao["tres"])*(admissao["pres"]/(admissao["epson"]*resultados["p1"] - admissao["pres"]))
    resultados["nres"] = comb_coef["reagentes"]*resultados["efe"]
    resultados["t1"]  = (admissao["Tamb"] + admissao["deltaT"] + resultados["efe"]*admissao["tres"])/(1 + resultados["efe"])
    resultados["rendimento_vol"] = (resultados["p1"]/admissao["pamb"])*(admissao["epson"]/(admissao["epson"]-1))*(admissao["Tamb"]/(resultados["t1"]*(1+resultados["efe"])))
    # %%        Passo 3: Compressao

    compressao = {"expo": compexo}        # Varia entre 1.3 e 1.37 

    resultados["t2"] = resultados["t1"]*(admissao["epson"]**(compressao["expo"] - 1))
    resultados["t2_c"] = resultados["t2"] - 273.15
    resultados["p2"] = resultados["p1"]*(admissao["epson"]**compressao["expo"])

    compressao["ei"] = {             # Energia interna
        "CO2": -0.43922+0.032813*(resultados["t2_c"])+0.00001036*(resultados["t2_c"]**2)-0.0000000018145*(resultados["t2_c"]**3),
        "CO": 0.00939989+0.0200386*(resultados["t2_c"])+0.00000389639*(resultados["t2_c"]**2)-0.000000000605288*(resultados["t2_c"]**3),
        "H2O": 0.146841+0.0231877*(resultados["t2_c"])+0.00000786405*(resultados["t2_c"]**2)-0.000000000829618*(resultados["t2_c"]**3),
        "H2": 0.101+0.01977865*(resultados["t2_c"])+0.000001602479*(resultados["t2_c"]**2)+0.0000000000429123*(resultados["t2_c"]**3),
        "N2": 0.03869+0.01969*(resultados["t2_c"])+0.000003708*(resultados["t2_c"]**2)-0.00000000055124*(resultados["t2_c"]**3),
        "O2": -0.133567+0.02181*(resultados["t2_c"])+0.00000355345*(resultados["t2_c"]**2)-0.00000000049*(resultados["t2_c"]**3),
    }


    ei = compressao["ei"]
    compressao["gases"] = comb_coef["CO2"] + comb_coef["CO"] + comb_coef["H2O"] + comb_coef["H2"] + comb_coef["N2"]
    compressao["ei"]["gases"] = (ei["CO2"]*comb_coef["CO2"] + ei["CO"]*comb_coef["CO"] + ei["H2O"]*comb_coef["H2O"] + ei["H2"]*comb_coef["H2"] + ei["N2"]*comb_coef["N2"])/compressao["gases"]
    compressao["ei"]["ar"] = (ei["O2"] + 3.76*ei["N2"])/4.76

    # %%        Passo 4: Combustao

    combustao = {
        "beta": 0.9,
        "pci": PCI[fuel["tipo"]],
        "u2": compressao["ei"]["ar"],
        "ulin2": compressao["ei"]["gases"],
        "delcl": ((1-fuel["alfa"])*fuel["peso"]*0.52193),
        "mi": (comb_coef["produtos"] + resultados["nres"])/(comb_coef["reagentes"] + resultados["nres"])
    }
    combustao["gases"] = compressao["gases"]

    param1 = (combustao["beta"]*(combustao["pci"] - combustao["delcl"]))/(comb_coef["reagentes"]*(1 + resultados["efe"]))
    param2 = (combustao["u2"] + (resultados["efe"]*combustao["ulin2"]))/(1 + resultados["efe"])
    u3 = (param1 + param2)/combustao["mi"]

    digitar_c = digitar_c        # Validar o que é esta merda
    resultados["t3"] = digitar_c + 273.15
    resultados["p3"] = (combustao["mi"]*resultados["p2"]*resultados["t3"])/resultados["t2"]
    resultados["pmax"] = 0.85*resultados["p3"]

    Lambda = resultados["p3"]/resultados["p2"]

    combustao["ei"] = {
        "CO2": -0.43922 + 0.032813*(digitar_c) + 0.00001036*(digitar_c**2) - 0.0000000018145*(digitar_c**3),
        "CO": 0.00939989 + 0.0200386*(digitar_c) + 0.00000389639*(digitar_c**2) - 0.000000000605288*(digitar_c**3),
        "H2O": 0.146841 + 0.0231877*(digitar_c) + 0.00000786405*(digitar_c**2) - 0.000000000829618*(digitar_c**3),
        "H2": 0.101 + 0.01977865*(digitar_c) + 0.000001602479*(digitar_c**2) + 0.0000000000429123*(digitar_c**3),
        "N2": 0.03869 + 0.01969*(digitar_c) + 0.000003708*(digitar_c**2) - 0.00000000055124*(digitar_c**3),
        "O2": -0.133567 + 0.02181*(digitar_c) + 0.00000355345*(digitar_c**2) - 0.00000000049*(digitar_c**3)
    }

    ei = combustao["ei"]
    combustao["ei"]["gases"] = (ei["CO2"]*comb_coef["CO2"] + ei["CO"]*comb_coef["CO"] + ei["H2O"]*comb_coef["H2O"] + ei["H2"]*comb_coef["H2"] + ei["N2"]*comb_coef["N2"])/combustao["gases"]
    combustao["ei"]["ar"] = (ei["O2"] + 3.76*ei["N2"])/4.76

    delta_u3 = u3 - combustao["ei"]["gases"]

    # %%        Passo 5: Expansao

    expansao = {"expo": expexpo}       # Varia de 1.23 a 1.30
    expansao["tx_compr"] = taxa_compressao

    resultados["p4"] = resultados["p3"]*(1/expansao["tx_compr"]**expansao["expo"])
    resultados["t4"] = resultados["t3"]*(1/expansao["tx_compr"]**(expansao["expo"]-1))

    # %%        Passo 6: Pressao media indicada [bar]

    passo6 = {
        "fii": fii,                # Varia de 0.92 a 0.97
        "fii_gas": fii_gas             # Varia de 0.75 a 0.90
    }

    fator_1 = resultados["p1"]*(admissao["epson"]**compressao["expo"])/(admissao["epson"]-1)
    fator_2 = (Lambda/(expansao["expo"]-1))*(1-(1/(admissao["epson"]**(expansao["expo"]-1))))
    fator_3 = (1/(compressao["expo"]-1))*(1-(1/admissao["epson"]**(compressao["expo"]-1)))

    resultados["pmi"] = fator_1*(fator_2 - fator_3)
    resultados["prim_cor"] = resultados["pmi"]*passo6["fii"]
    resultados["delpgas"] = passo6["fii_gas"]*(admissao["pres"]-resultados["p1"])
    resultados["seg_cor"] = resultados["prim_cor"] - resultados["delpgas"]

    # %%       Passo 7: Pressao media de atrito

    arques = (0.5 + 0.18*(rotacao/1000) + 0.02*((rotacao/1000)**2))

    aux1 = 16.3761 + 2.28629*((rotacao/1000)) + 0.297053*(rotacao/1000)**2
    aux2 = 0.01*(5.44659 - 0.02495*(rotacao/1000) - 0.174376*(rotacao/1000)**2)
    abnt = (6.89*aux1 - aux2*resultados["seg_cor"])/((1 - aux2)*100)

    aux1 = 0.05
    aux2 = 0.0155
    khovakh = (aux1+(aux2*2*rotacao*curso/60))*10

    resultados["pma"] = {
        "arques": arques,
        "abnt": abnt,
        "khovakh": khovakh
    }

    # %%        Passo 8: Potencia efetiva

    Wi = (cilindrada*resultados["seg_cor"]*rotacao)/(0.6*x)
    arques = (resultados["seg_cor"]-resultados["pma"]["arques"])*cilindrada*rotacao/(0.6*x)
    abnt = (resultados["seg_cor"]-resultados["pma"]["abnt"])*cilindrada*rotacao/(0.6*x)
    khovakh = (resultados["seg_cor"]-resultados["pma"]["khovakh"])*cilindrada*rotacao/(0.6*x)

    resultados["potencia_efetiva"] = {
        "arques": arques,
        "abnt": abnt,
        "khovakh": khovakh
    }

    # Calculados a partir do ajuste derivado de Arques
    resultados["torque"] = 1000*100*(resultados["seg_cor"]-resultados["pma"]["arques"])*(diametro**2)/8*curso*num_cilindros/x
    resultados["rendimento_mec"] = (resultados["seg_cor"]-resultados["pma"]["arques"])/resultados["seg_cor"]

    # %%        Passo 9: Trabalho indicado por ciclo, por kg de combustivel

    R = 8.314
    V1 = comb_coef["reagentes"]*R*resultados["t1"]/(resultados["p1"]*100)
    V2 = V1/taxa_compressao
    pmi = resultados["seg_cor"]*100

    resultados["Wi"] = pmi*(V1-V2)/1000
    resultados["rendimento_ind"] = resultados["Wi"]/combustao["pci"]
    resultados["rendimento_ter"] = resultados["rendimento_vol"] * resultados["rendimento_mec"] * resultados["rendimento_ind"]
    resultados["cc"] = resultados["potencia_efetiva"]["arques"]/combustao["pci"]/resultados["rendimento_ter"]
    resultados["cec"] = resultados["cc"]*3600/resultados["potencia_efetiva"]["arques"]

    #print(resultados)
    # %%        Resultados Finais
    ref_pe = referencia["potencia_efetiva"]
    res_pe = resultados["potencia_efetiva"]["arques"]/0.7351    # Em HP
    ref_tq = referencia["torque"]
    res_tq = resultados["torque"]

    return res_pe, res_tq, resultados["rendimento_ter"], resultados["cec"]

def paraPorcentagem(lista):
    y0 = lista[0]
    for i in range(len(lista)):
        lista[i] = (lista[i]-y0)*100/y0
    return lista

# Variação de temperatura e pressão atmosférica para potencia

range_tamb = np.arange(250, 350, 0.1)
potencia, torque, rendTermico, cec = [], [], [], []
for tamb in range_tamb:
    p, t, r, c = simulacao(tamb, 1.013, op="potencia", printTable=0)
    potencia.append(p); torque.append(t); rendTermico.append(r); cec.append(c)

fig, ax = plt.subplots(figsize=(8, 5))
ln1 = ax.plot(range_tamb, paraPorcentagem(potencia), label="Potência Efetiva")
ln2 = ax.plot(range_tamb, paraPorcentagem(rendTermico), label="Rendimento Térmico")
ln3 = ax.plot(range_tamb, paraPorcentagem(cec), label="CEC")
lns = ln1 + ln2 + ln3
labs = [l.get_label() for l in lns]
plt.legend(labs)
plt.grid(linestyle="--")
plt.ylabel("Variação [%]")
plt.xlabel("Temperatura ambiente [°C]")
plt.show()

range_pamb = np.arange(900, 1100, 0.1)
potencia, torque, rendTermico, cec = [], [], [], []
for pamb in range_pamb:
    p, t, r, c = simulacao(298, pamb, op="potencia", printTable=0)
    potencia.append(p); torque.append(t); rendTermico.append(r); cec.append(c)

fig, ax = plt.subplots(figsize=(8, 5))
ln1 = ax.plot(range_pamb, paraPorcentagem(potencia), label="Potência Efetiva")
ln2 = ax.plot(range_pamb, paraPorcentagem(rendTermico), label="Rendimento Térmico")
ln3 = ax.plot(range_pamb, paraPorcentagem(cec), label="CEC")
lns = ln1 + ln2 + ln3
labs = [l.get_label() for l in lns]
plt.legend(labs)
plt.grid(linestyle="--")
plt.ylabel("Variação [%]")
plt.xlabel("Pressão Atmosférica [kPa]")
plt.show()

# Variação de temperatura e pressão atmosférica para torque

range_tamb = np.arange(250, 350, 0.1)
potencia, torque, rendTermico, cec = [], [], [], []
for tamb in range_tamb:
    p, t, r, c = simulacao(tamb, 1.013, op="torque", printTable=0)
    potencia.append(p); torque.append(t); rendTermico.append(r); cec.append(c)

fig, ax = plt.subplots(figsize=(8, 5))
ln1 = ax.plot(range_tamb, paraPorcentagem(torque), label="Torque")
ln2 = ax.plot(range_tamb, paraPorcentagem(rendTermico), label="Rendimento Térmico")
ln3 = ax.plot(range_tamb, paraPorcentagem(cec), label="CEC")
lns = ln1 + ln2 + ln3
labs = [l.get_label() for l in lns]
plt.legend(labs)
plt.grid(linestyle="--")
plt.ylabel("Variação [%]")
plt.xlabel("Temperatura ambiente [°C]")
plt.show()

range_pamb = np.arange(900, 1100, 0.1)
potencia, torque, rendTermico, cec = [], [], [], []
for pamb in range_pamb:
    p, t, r, c = simulacao(298, pamb, op="torque", printTable=0)
    potencia.append(p); torque.append(t); rendTermico.append(r); cec.append(c)

fig, ax = plt.subplots(figsize=(8, 5))
ln1 = ax.plot(range_pamb, paraPorcentagem(torque), label="Torque")
ln2 = ax.plot(range_pamb, paraPorcentagem(rendTermico), label="Rendimento Térmico")
ln3 = ax.plot(range_pamb, paraPorcentagem(cec), label="CEC")
lns = ln1 + ln2 + ln3
labs = [l.get_label() for l in lns]
plt.legend(labs)
plt.grid(linestyle="--")
plt.ylabel("Variação [%]")
plt.xlabel("Pressão Atmosférica [kPa]")
plt.show()