import torch
import numpy as np
import matplotlib.pyplot as plt
from torch.autograd import Variable

import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.model_selection import train_test_split
# %matplotlib inline

import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import math
import torchvision
import torchvision.transforms as transforms
from torchvision.utils import make_grid
import scipy.io as spio
from os.path import dirname, join as pjoin
import os,sys
import torch.optim.lr_scheduler as lrs
from torch.optim import lr_scheduler


# data = [[1,2],[3,4]]
# tensor = torch.FloatTensor(data) #32-bit floating point 

# print(
#     '\nnumpy:', np.matmul(data, data),
#     '\ntorch:', torch.mm(tensor, tensor)
# )
# data = np.array(data)
# print(
#     '\nnumpy:',data.dot(data),
#     '\ntorch:',tensor.dot(tensor)
# )

# tensor = torch.FloatTensor([[1,2], [3,4]])
# variable = Variable(tensor, requires_grad = True)
# print(tensor)
# print(variable)

# t_out = torch.mean(tensor*tensor)
# v_out = torch.mean(variable*variable)
# print(t_out)
# print(v_out)

# v_out.backward()
# print(variable.grad)
# print(variable.data)

# x = torch.linspace(-5, 5, 200)
# x = Variable(x)
# x_np = x.data.numpy()

# y_relu = torch.relu(x).data.numpy()
# y_sigmoid = torch.sigmoid(x).data.numpy()
# y_tanh = torch.tanh(x).data.numpy()
# y_softplus = torch.softplus(x).data.numpy()
# y_softmax = torch.softmax(x)

# plt.figure(1, figsize = (8, 6))
# plt.subplot(221)
# plt.plot(x_np, y_relu, c='red', label='relu')
# plt.ylim((-1, 5))
# plt.legend(loc = 'best')

# plt.subplot(222)
# plt.plot(x_np, y_sigmoid, c='red', label='sigmoid')
# plt.ylim((-0.2, 1.2))
# plt.legend(loc = 'best')

# plt.subplot(223)
# plt.plot(x_np, y_tanh, c='red', label='tanh')
# plt.ylim((-1.2, 1.2))
# plt.legend(loc = 'best')
# plt.show()

# plt.subplot(224)
# plt.plot(x_np, y_softplus, c='red', label='softplus')
# plt.ylim((-1.2, 1.2))
# plt.legend(loc = 'best')

# x = torch.unsqueeze(torch.linspace(-1, 1, 100), dim = 1) # x data(tensor), shape=(100, 1)
# y = x.pow(2) + 0.2*torch.rand(x.size())

filename=os.path.basename(os.path.realpath(sys.argv[0]))
dirname=os.path.dirname(os.path.realpath(sys.argv[0]))

M, input_size, hidden_size, output_size = 400, 2, 100, 2
data_dir = pjoin(dirname)
# mat_fname = pjoin(data_dir, 'rand_10000.mat')
mat_fname = pjoin(data_dir, 'cpp_random_sample_50000.mat')

mat_data = spio.loadmat(mat_fname)
q = mat_data['q']
x_train = torch.from_numpy(q).float()
# scale units
# x_train_max, _ = torch.max(x_train, 0)
# x_train = torch.div(x_train, x_train_max)


Im = mat_data['JMJ']
Im = Im[:,2:3]
y_train = torch.from_numpy(Im).float()
# scale units
y_train_max, _ = torch.max(y_train, 0)
y_train = torch.div(y_train, y_train_max)

# mat_fname = pjoin(data_dir, 'rand_900.mat')
mat_fname = pjoin(data_dir, 'cpp_test_4500.mat')

mat_data = spio.loadmat(mat_fname)
q = mat_data['q']
x_test = torch.from_numpy(q).float()
# x_test_max, _ = torch.max(x_test, 0)
# x_test = torch.div(x_test, x_test_max)

Im = mat_data['JMJ']
Im = Im[:,2:3]
y_test = torch.from_numpy(Im).float()
# y_test_max, _ = torch.max(y_test, 0)
y_test = torch.div(y_test, y_train_max)


if torch.cuda.is_available():
    x_train = x_train.cuda()
    x_test = x_test.cuda() 
    y_train = y_train.cuda()
    y_test = y_test.cuda() 

def euclidean(x0, x1):
    x0, x1 = np.array(x0), np.array(x1)
    d = np.sum((x0 - x1)**2)**0.5
    return d

# class Net(torch.nn.Module):
#     def __init__(self, n_feature, n_hidden, n_output):
#         super(Net, self).__init__()
#         self.hidden = torch.nn.Linear(n_feature, n_hidden)
#         self.predict = torch.nn.Linear(n_hidden, n_output)
    
#     def forward(self, x):
#         x = torch.relu(self.hidden(x))
#         x = self.predict(x)
#         return x

# net = Net(2, 100, 2)
# print(net)

# plt.ion() 
# plt.show()

# optimizer = torch.optim.SGD(net.parameters(), lr=0.005)
# loss_func = torch.nn.MSELoss()
# plt.scatter(x.data.numpy(), y.data.numpy())

# for t in range(0):
#     prediction = net(x)
#     loss = loss_func(prediction, y)

#     optimizer.zero_grad()
#     loss.backward()
#     optimizer.step()
#     if t % 5 == 0:
#         # plot and show learning process
#         plt.cla()
#         plt.scatter(x.data.numpy(), y.data.numpy())
#         plt.plot(x.data.numpy(), prediction.data.numpy(), 'r-', lw = 5)
#         plt.text(0.5, 0, 'Loss=%.4f' % loss.data.item(), fontdict={'size':20, 'color': 'red'})
#         plt.pause(0.1)

# plt.ioff()
# plt.show()


# net2 = torch.nn.Sequential(
#     torch.nn.Linear(2,100),
#     torch.nn.ReLU(),
#     torch.nn.Linear(100,2)
# )

# print(net)
# print(net2)

def save():
    # save net1
    input_size = 2
    hidden_size = 100
    output_size = 2
    learning_rate = 0.002
    momentum_ = 0.9
    epoch_ = 50000
    # Use torch.jit.trace to generate a torch.jit.ScriptModule via tracing.

    net1 = torch.nn.Sequential(
        torch.nn.Linear(input_size, hidden_size),
        torch.nn.Tanh(),
        torch.nn.Linear(hidden_size, hidden_size), #1
        torch.nn.Tanh(),
        torch.nn.Linear(hidden_size, hidden_size), #2
        torch.nn.ReLU(),
        torch.nn.Linear(hidden_size, hidden_size), #3
        torch.nn.ReLU(),
        torch.nn.Linear(hidden_size, hidden_size), #4
        torch.nn.ReLU(),
        torch.nn.Linear(hidden_size, hidden_size), #5
        torch.nn.ReLU(),
        # torch.nn.Linear(hidden_size, hidden_size), #6
        # torch.nn.ReLU(),
        # torch.nn.Linear(hidden_size, hidden_size), #7
        # torch.nn.ReLU(),
        # torch.nn.Linear(hidden_size, hidden_size), #8
        # torch.nn.ReLU(),
        # torch.nn.Linear(hidden_size, hidden_size), #9
        # torch.nn.ReLU(),
        # torch.nn.Linear(hidden_size, hidden_size), #10
        # torch.nn.ReLU(),
        # torch.nn.Linear(hidden_size, hidden_size), #11
        # torch.nn.ReLU(),
        torch.nn.Linear(hidden_size, output_size)
    )
    
    if torch.cuda.is_available():
        net1 = net1.cuda()

    # optimizer = torch.optim.Adam(net1.parameters(), lr = learning_rate)
    optimizer = torch.optim.SGD(net1.parameters(), lr = learning_rate, momentum = momentum_)
    # Decay LR by a factor of 0.1 every 100 epochs
    # scheduler = lr_scheduler.StepLR(optimizer, step_size=epoch_/3, gamma=0.1)
    loss_func = torch.nn.MSELoss()
    for t in range(epoch_):      
        prediction = net1(x_train)
        loss = loss_func(prediction, y_train)

        prediction_test = net1(x_test)
        accuracy = loss_func(prediction_test, y_test)

        if t % 100 == 0:
            print('Epoch {}/{}'.format(t, epoch_ - 1))
            print('-' * 10)
            print('loss', loss)
            print('acc', accuracy)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        # scheduler.step()

    save_path = 'D:/Research/Muscles/Projects/DeepLearning'
    torch.save(net1, os.path.join(save_path, "net.pkl"))    
    torch.save(net1, 'net_.pkl') #entire net
    # torch.save(net1.state_dict(), 'net_params.pkl') # parameters


def restore_net():
    save_path = 'D:/Research/Muscles/Projects/DeepLearning'
    name_of_file = 'trained_900'
    completeName = os.path.join(save_path, name_of_file+".pt")         


    net2 = torch.load(os.path.join(save_path, "net.pkl"))
    if torch.cuda.is_available():
        example = torch.rand(1, 2).cuda()
    else:
        example = torch.rand(1, 2)    
    # 
    

    traced_script_module = torch.jit.trace(net2, example)
    traced_script_module.save(completeName)
    # output = traced_script_module(torch.ones(1,2))
    tensor = torch.FloatTensor([[0.2609, -1.5412]])
    if torch.cuda.is_available():
        output = traced_script_module(tensor.cuda())
    else:
        output = traced_script_module(tensor)
        print("cpu")
    # output = traced_script_module(torch.ones(1,2).cuda())
    print(output)

def restore_params():
    net3 = torch.load('net.pkl')


    # net3.load_state_dict(torch.load('net_params.pkl'))
    optimizer = torch.optim.SGD(net3.parameters(), lr = 0.01)
    loss_func = torch.nn.MSELoss()

    for t in range(50000):
        prediction = net3(x_train)
        loss = loss_func(prediction, y_train)

        prediction_test = net3(x_test)
        accuracy = loss_func(prediction_test, y_test)

        if t % 20 == 0:
            print('loss', loss)
            print('acc', accuracy)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    
    torch.save(net3, 'net.pkl') #entire net
    torch.save(net3.state_dict(), 'net_params.pkl') # parameters

    # prediction = net3(x_test)
    # loss = loss_func(prediction, y_test)
    # print('test')
    # print(loss)

save()
#restore_net()
#restore_params()