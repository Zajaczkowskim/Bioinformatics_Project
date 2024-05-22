from model import Model

model = Model(500,20)
model.run_simulation()

print("Expected loops: ")
print(model.expected_loops)