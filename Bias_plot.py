import numpy as np
import matplotlib.pyplot as plt

def main():
    year=['16','17','18']
    mass=['M5','M15','M30']

    bias_test_16M5=[0.1156,-0.05385,0.02754]#Bias for each test function
    bias_name_16M5=['5th Order Bernstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_16M5=[9.277]

    bias_test_16M15=[-0.3492,-0.03918,-0.0822]#Bias for each test function
    bias_name_16M15=['5th Order Bernstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_16M15=[12.646]

    bias_test_16M30=[0.1747,-0.005656,-0.01261]#Bias for each test function
    bias_name_16M30=['5th Order Bernstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_16M30=[12.012]
    #
    bias_test_17M5=[0.1097,0.03437,-0.05286,-0.009763]#Bias for each test function
    bias_name_17M5=['5th Order Bernstein', '6th Order Berstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_17M5=[5.078]

    bias_test_17M15=[-0.2327,-0.1479,0.02248,-0.04551]#Bias for each test function
    bias_name_17M15=['5th Order Bernstein', '6th Order Berstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_17M15=[7.471]

    bias_test_17M30=[0.1297,0.1274,-0.01048,-0.002183]#Bias for each test function
    bias_name_17M30=['5th Order Bernstein', '6th Order Berstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_17M30=[6.689]
    #
    bias_test_18M5=[0.09222,-0.1148,-0.03817]#Bias for each test function
    bias_name_18M5=['5th Order Bernstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_18M5=[5.029]

    bias_test_18M15=[-0.2593,0.02349,-0.07739]#Bias for each test function
    bias_name_18M15=['5th Order Bernstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_18M15=[6.567]

    bias_test_18M30=[0.1851,-0.03245,-0.04167]#Bias for each test function
    bias_name_18M30=['5th Order Bernstein', '3rd Order Exponential', '3rd Order Power Low']
    bias_true_18M30=[5.859]

    bias_test_all=[bias_test_16M5, bias_test_16M15, bias_test_16M30, bias_test_17M5, bias_test_17M15, bias_test_17M30, bias_test_18M5, bias_test_18M15, bias_test_18M30]
    bias_name_all=[bias_name_16M5, bias_name_16M15, bias_name_16M30, bias_name_17M5, bias_name_17M15, bias_name_17M30, bias_name_18M5, bias_name_18M15, bias_name_18M30]
    bias_true_all=[bias_true_16M5, bias_true_16M15, bias_true_16M30, bias_true_17M5, bias_true_17M15, bias_true_17M30, bias_true_18M5, bias_true_18M15, bias_true_18M30]

    for i in range(len(bias_test_all)):
        bias_test=bias_test_all[i]
        bias_name=bias_name_all[i]
        bias_true=bias_true_all[i]

        x=np.array(bias_test)
        y=0.5*np.ones(len(bias_test))

        plt.xlabel('$\\frac{\mu-\\tilde{\mu}}{\sigma_{\mu}}$')
        plt.ylabel('$\\tilde{\mu}$')
        plt.xlim(xmax=0.5,xmin=-0.5)
        plt.ylim(ymax=1.0,ymin=0)

        colors = ['red','orange','blue','purple','darkturquoise']
        area = np.pi * 4**2

        plt.plot([0.14, 0.14], [0., 1], c='black', linestyle='--')
        plt.plot([-0.14, -0.14], [0., 1], c='black', linestyle='--')
        plt.yticks(np.array([0.5]),np.array(bias_true))
        for j in range(len(bias_test)):
            plt.scatter(x[j], y[j], s=area, c=colors[j], alpha=0.4, label=bias_name[j])

        plt.title('Toy Function: 3rd Order Power Low (' + year[int(i/3)] + '_'+mass[i%3] + ')')
        plt.grid()
        plt.legend()
        plt.savefig('./BiasPlot/' + year[int(i/3)] + '_'+mass[i%3] + '.png')
        plt.close('all')
        #plt.show()

main()
