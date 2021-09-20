# This is a sample Python script.


import matplotlib.pyplot as plt
import numpy as np
import os
import PIL
import pathlib

from tensorflow.keras import datasets, layers, models
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential


# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')


data_dir = pathlib.Path("/home/kostasl/workspace/zebrafishtrack/img/trainset/")

fish = list(data_dir.glob('./fish/*'))
PIL.Image.open(str(fish[0]))

nonfish = list(data_dir.glob('./nonfish/*'))
PIL.Image.open(str(nonfish[0]))

batch_size = 32
img_height = 28
img_width = 38

train_ds = tf.keras.preprocessing.image_dataset_from_directory(
  data_dir,
  validation_split=0.2,
  subset="training",
  seed=123,
  image_size=(img_height, img_width),
  batch_size=batch_size)



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
