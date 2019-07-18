import tensorflow as tf
from tensorflow import keras

from keras.models import Sequential
from keras.layers import Dense, Activation

model = Sequential()
# Adds to the model a densely-connected layer with 27 units with input shape 2, an (x,y) pair:
model.add(Dense(27, input_shape=(2,), activation='relu'))
# Adds another layer with 18 units, each connected to 27 outputs of previous layer
model.add(Dense(18, activation='relu'))
# Last layer with 9 units, each connected to 18 outputs of previous layer
model.add(Dense(9, activation='softmax'))
