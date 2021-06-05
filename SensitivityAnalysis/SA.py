

def SA_morris(sample=False,analysis=False):
    from SALib.sample import morris as morrisS
    from SALib.analyze import morris
    from SALib.test_functions import Ishigami
    import numpy as np

    # Define the model inputs
    problem = {
    'num_vars': 2, # Number of parameters
    'names': ['par1', 'par2'], # Name of parameters
    'bounds': [[0, 1], # Bounds of parameters
    [50, 100]]
    }
    
    if (sample == True):
        # Generate samples
        param_values = morrisS.sample(problem, 6)
        print(param_values.shape)
        np.save('SA_samples', param_values)        
    
    if (analysis == True):
        param_values = np.load('SA_samples.npy')
        # Run model (example)
        Ntimes = 10
        Y = np.zeros((param_values.shape[0], Ntimes))    
        for i, param_value in enumerate(param_values):
            Y[i,:] = Model(param_value,Ntimes)
        print(Y.shape)

        # Perform analysis
        Si = []
        for k in range(Ntimes):
            Si.append(morris.analyze(problem, param_values, Y[:,k], conf_level=0.95, print_to_console=False))

        # Print the first-order sensitivity indices
        for IndexSens in Si:
            print('first-order sensitivity Mu_star: ' + str(IndexSens['mu_star']))

if __name__ == '__main__':
    print("\n\n Sensitivity analysis - MORRIS METHOD\n")
    SA_morris(sample=True)
