import re
import numpy as np
import glob
import matplotlib.pyplot as plt

def extract_number(file_name):
    match = re.search(r"velocity_([\d\.]+)_", file_name)
    return float(match.group(1)) if match else float('inf')  # Use a high number for unmatched cases


class Data:
    """
    Data for reading and plotting from OF-TCAT
    """
    def __init__(self,files, num_variables):
        self.files = files
        self.num_files = len(files)
        self.num_variables = num_variables
        self.data = self.initialize()

    def initialize(self):
        """
            Initialize the numpy arrays for each variable
        """

        data = {}
        with open(self.files[0], "r", encoding="utf-8") as f:
            lines = f.readlines()
            for n in range(self.num_variables):
                split = lines[n].split(",")
                name = split[0]
                if len(split) == 2:
                    value = split[-1].split("\n")[0]
                    units = None
                else:
                    value = split[-2].split("\n")[0]
                    units = split[-1].split("\n")[0]
                if "(" in value:
                    num_var = len(value.split(" "))
                    data[name] = {'dim':num_var,'data':np.zeros([self.num_files,num_var]),'units':units}
                else:
                    data[name] = {'dim':1,'data':np.zeros(self.num_files),'units':units}

        return data

    def get_file_id(self,file):
        """
        Determine the id of a file 

        Args:
            file (_type_): _description_
        """
        f_split = file.split("_")[-2]
        n = f_split.split(".txt")[0]

        return int(n) 

    def read_output(self):
        """
        Read the output from OF-TCAT

        Assume not sorted!!!

        Args:
            file (_type_): _description_
        """
        for n_file,file in enumerate(self.files):
            with open(file, "r", encoding="utf-8") as f:
                lines = f.readlines()
                print(n_file,file)
                for n in range(self.num_variables):
                    split = lines[n].split(",")
                    name = split[0]
                    
                    if len(split) == 2:
                        value = split[-1].split("\n")[0]
                    else:
                        value = split[-2].split("\n")[0]
                    
                    if self.data[name]['dim'] > 1:
                        numbers = value.strip('( )').split()
                        out_value = np.array([float(num) for num in numbers])
                        self.data[name]['data'][n_file] = out_value
                    else:                
                        self.data[name]['data'][n_file] = value



folder_in = "tcat/"
files = glob.glob(folder_in+"*")
sorted_files = sorted(files, key=extract_number)

steady_data =  Data(sorted_files,28)
steady_data.read_output()

Re = steady_data.data['Reynolds Number']['data'][:]
rho = steady_data.data['Macroscale rho']['data'][:]
U = steady_data.data['Macroscale U']['data'][:,0]

mu_w = steady_data.data['Macroscale e_w*rho*grad(chem potential)']['data'][:,0]
psi_w = steady_data.data['Macroscale ew*rho*grad(potential)']['data'][:,0]

### Error Plots

plt.loglog(Re,np.abs(steady_data.data['Microscale mass error']['data'][:]),'o', label = "Mass")
plt.loglog(Re,np.abs(steady_data.data['Microscale momentum error']['data'][:,0]),'o',label = "Momentum")
plt.loglog(Re,np.abs(steady_data.data['Macroscale momentum error']['data'][:,0]),'o',label = "Momentum")
plt.xlabel("Reynolds Number")
plt.ylabel("Absolute Residual")
plt.legend()
plt.show()

# ### Term Plots
plt.semilogx(Re,steady_data.data['Macroscale ddt(rho*U)']['data'][:,0],'o',label = 'ddt(rho*U)')
plt.semilogx(Re,steady_data.data['Macroscale div(rho*U*U)']['data'][:,0],'o',label = 'div(rho*U*U)')
plt.semilogx(Re,steady_data.data['Macroscale rho*g']['data'][:,0],'o',label = 'rho*g')
plt.semilogx(Re,steady_data.data['Macroscale div(t)']['data'][:,0],'o',label = 'div(t)')
plt.semilogx(Re,steady_data.data['Macroscale T']['data'][:,0],'o',label = 'T')
plt.xlabel("Reynolds Number")
plt.ylabel("Macroscale Momentum Term")
plt.legend()
plt.show()

### Potential Gradient Plots

plt.loglog(Re,mu_w,'o', label = "grad mu")
plt.loglog(Re,psi_w,'o', label = "grad psi")
plt.xlabel("Reynolds Number")
plt.ylabel("U")
plt.legend()
plt.show()


### Darcy
plt.loglog(U/rho,-(mu_w+psi_w)/(rho*U*U),'o')
plt.xlabel("Reynolds Number")
plt.ylabel("Potential")
plt.show()

# ### Resistance Tensor
# plt.loglog(Re,np.abs(steady_data.data['K']['data'][:]),'o')
# plt.xlabel("Reynolds Number")
# plt.ylabel("Resistance Tensor xx")
# plt.show()