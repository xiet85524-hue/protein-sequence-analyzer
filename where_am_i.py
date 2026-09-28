import os
print("script file :", __file__)
print("i am standing in :", os.getcwd())
print("output goes to :", os.path.dirname(os.path.abspath(__file__)))