from model import Model

model = Model(500,2)
model.run_simulation()
 
print("Expected loops: ")
print(model.expected_loops)