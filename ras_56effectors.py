from scipy.optimize import fsolve
from numpy import *
import pandas as pd

# Run function main().
# NB. Modify the file paths in main() to your local folder e.g. C:/Users/...

def system(p):
    """
    System of steady state equations.
    :return: list of the equations
    """

    C1, C2, C3, C4, C5, C6, C7, C8, \
    C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22, \
    C23, C24, C25, C26, C27, C28, C29, C30, C31, C32, C33, C34,\
    C35, C36, C37, C38, C39, C40, C41, C42, C43, \
    C44, C45, C46, C47, C48, C49, C50, C51, C52, C53, C54, C55, C56 = p

    C = [C1, C2, C3, C4, C5, C6, C7, C8,
         C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22,
         C23, C24, C25, C26, C27, C28, C29, C30, C31, C32, C33, C34, C35, C36, C37,
         C38, C39, C40, C41, C42, C43,
         C44, C45, C46, C47, C48, C49, C50, C51, C52, C53, C54, C55, C56]

    eqs = [C[i] * (Kd[i] + Rtot - sum(C)) + Etot[i] * (sum(C) - Rtot) for i in range(n)]

    return eqs


def params_and_sols(Kd_microM, Rtot, Etot):
    """
    Solve the system of equations for given parameters.
    """

    global n, Kd
    n = len(Kd_microM)  # nbr of effectors

    # micro to nano Molar transformation
    micro_to_nano = [10 ** 3] * n
    Kd = multiply(Kd_microM, micro_to_nano)

    # Assign ICs to solve the system
    IC = [1.] * n

    # Compute absolute values for complex formation (nM)
    C = fsolve(system, IC)
    # Compute relative/normalized values for complex formation (%)
    C_perc = [C[i] / sum(C) for i in range(n)]

    # verify if the results are 0:
    if [C[i] * (Kd[i] + Rtot - sum(C)) + Etot[i] * (sum(C) - Rtot) for i in range(n)] > [10 ** (-2)] * n:
        print('ATTENTION system not solved properly! \n')
        # break

    return C_perc, C


def main(GTP=.2, MUT=False):
    """
    Compute Ras-effector complex concentrations based on the GTP load of H,K,NRAS.
    Save the data into .xlsx files.
    Usually set GTP=.2 if  MUT is True.

    :param GTP: .2, .5, .75, .9, 1
    :param MUT: True, False
    :return: dataframes of the computed complex concentrations in percent and nanoMolar
    """

    global Rtot, Etot
    if MUT:
        print('Choose one mutant (type the associated number): HRAS (0), KRAS (1), NRAS (2)')
        WHICH_MUT = input('> ')
        WHICH_MUT = int(WHICH_MUT)
        string_mut = ['H', 'K', 'N']
        ras_mut = string_mut[WHICH_MUT]

        print(f'Choose GTP load for {ras_mut}RAS mutant: e.g. .5, .75, 1, 1.25, 1.5')
        MUT_GTP = input('> ')
        MUT_GTP = float(MUT_GTP)
    else:
        MUT_GTP = .2

    ras = [GTP]*3
    if MUT:
        ras[WHICH_MUT] = MUT_GTP

    print(f'GTP load of HRAS, KRAS, NRAS: {ras}')

    f = 'C:/Users/simo_/input_data_SCatozzi.xlsx'
    wb = pd.ExcelFile(f)  # load workbook
    ws = wb.sheet_names  # worksheet list
    dg = wb.parse(sheet_name=ws[0], header=0, index_col=None, usecols=None)
    dg.fillna(0, inplace=True)
    dg.rename(columns={'Kd (uM)': 'Kd'}, inplace=True)

    efflist = dg.Gene_symbol.tolist()[3:]  # effector list
    Kdlist = dg.Kd.tolist()[3:]  # omit first 3 N/A values for HKNRAS
    tissuelist = dg.columns[3:]

    # Dataframe with effectors data
    de = dg.loc[3:]
    # Dataframe with pan-RAS data: sum of all HKN-RAS
    dr = multiply(dg.iloc[0, 3:], ras[0]) + \
         multiply(dg.iloc[1, 3:], ras[1]) + \
         multiply(dg.iloc[2, 3:], ras[2])

    # Initialize dataframes
    df_perc = pd.DataFrame(index=efflist)
    df_perc['Class'] = de.Class.tolist()
    df_nM = pd.DataFrame(index=efflist)
    df_nM['Class'] = de.Class.tolist()

    # Compute percent and nanoMolar complex concentrations
    for t in tissuelist:
        Etot = de[t].tolist()
        Rtot = dr[t].tolist()

        C_perc, C_nM = params_and_sols(Kdlist, Rtot, Etot)

        # Fill in dataframes with percent or nanoMolar values
        df_perc[t] = multiply(C_perc, 100)
        df_nM[t] = C_nM

    # Save the dataframes as .xlsx files
    if not MUT:
        file_perc = f'C:/Users/simo_/Complexes_{int(GTP*100)}PanRAS_perc.xlsx'
        file_nM = f'C:/Users/simo_/Complexes_{int(GTP*100)}PanRAS_nM.xlsx'
    else:
        file_perc = f'C:/Users/simo_/Complexes_{int(MUT_GTP*100)}{ras_mut}RAS_MUT_perc.xlsx'
        file_nM = f'C:/Users/simo_/Complexes_{int(MUT_GTP*100)}{ras_mut}RAS_MUT_nM.xlsx'

    df_perc.to_excel(file_perc, sheet_name='C% 56 effectors')
    df_nM.to_excel(file_nM, sheet_name='C_nM 56 effectors')
    print(file_perc, 'saved.')
    print(file_nM, 'saved.')

    return df_perc, df_nM
