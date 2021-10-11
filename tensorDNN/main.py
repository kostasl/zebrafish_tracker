# This is a sample Python script.
import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np
import os
import io
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


print(tf.version.VERSION)


data_dir = pathlib.Path("/home/kostasl/workspace/zebrafishtrack/tensorDNN/trainset/")
valid_dir = pathlib.Path("/home/kostasl/workspace/zebrafishtrack/tensorDNN/trainset/")

fish = list(data_dir.glob('./fish/*.jpg'))
#PIL.Image.open(str(fish[0]))


nonfish = list(data_dir.glob('./nonfish/*.jpg'))
#PIL.Image.open(str(nonfish[0]))

batch_size = 64
img_height = 38
img_width = 28
def train_model(batch_size,img_height,img_width):

    train_ds = tf.keras.preprocessing.image_dataset_from_directory(
      str(data_dir),
      validation_split=0.01,
      subset="training",
      seed=123,
      image_size=(img_height, img_width),
      color_mode='grayscale',
      batch_size=batch_size)


    val_ds = tf.keras.preprocessing.image_dataset_from_directory(
      str(valid_dir),
      validation_split=0.9,
      subset="validation",
      seed=123,
      image_size=(img_height, img_width),
      color_mode='grayscale',
      batch_size=batch_size)

    class_names = train_ds.class_names
    print(class_names)


    plt.figure(figsize=(10, 10))
    for images, labels in train_ds.take(1):
      for i in range(9):
        ax = plt.subplot(3, 3, i + 1)
        plt.imshow(images[i].numpy().astype("uint8"))
        plt.title(class_names[labels[i]])
        plt.axis("off")

    for image_batch, labels_batch in train_ds:
      print(image_batch.shape)
      print(labels_batch.shape)
      break

    ## Data Augmentation
    data_augmentation = keras.Sequential(
      [
        layers.experimental.preprocessing.RandomFlip("horizontal",
                          input_shape=(img_height,  img_width, 1)),
        layers.experimental.preprocessing.RandomRotation(0.05),
        layers.experimental.preprocessing.RandomZoom(0.1),
      ]
    )

    plt.figure(figsize=(10, 10))
    for images, _ in train_ds.take(1):
      for i in range(9):
        augmented_images = data_augmentation(images)
        ax = plt.subplot(3, 3, i + 1)
        plt.imshow(augmented_images[0].numpy().astype("uint8"))
        plt.axis("off")


    ## Config Dataset performance
    AUTOTUNE = tf.data.AUTOTUNE

    train_ds = train_ds.cache().shuffle(1000).prefetch(buffer_size=AUTOTUNE)
    val_ds = val_ds.cache().prefetch(buffer_size=AUTOTUNE)

    ##Standardize the data
    normalization_layer = tf.keras.layers.experimental.preprocessing.Rescaling(1./255)

    ## Create the model ##
    num_classes = 2

    model = Sequential([
      data_augmentation,
      layers.experimental.preprocessing.Rescaling(1./255, input_shape=(img_height, img_width, 1)),
      layers.Conv2D(16, 3, padding='same', activation='relu'),
      layers.MaxPooling2D(),
      layers.Conv2D(32, 3, padding='same', activation='relu'),
      layers.MaxPooling2D(),
      layers.Conv2D(64, 3, padding='same', activation='relu'),
      layers.MaxPooling2D(),
      layers.Flatten(),
      layers.Dense(128, activation='relu'),
      layers.Dense(num_classes)
    ])
    ##COMPILE MODEL
    model.compile(optimizer='adam',
                  loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
                  metrics=['accuracy'])


    model.summary()

    checkpoint_path = "training_1/cp.ckpt"
    checkpoint_dir = os.path.dirname(checkpoint_path)

    # Create a callback that saves the model's weights
    cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath=checkpoint_path,
                                                     save_weights_only=True,
                                                     verbose=1)


    ## TRAIN MODEL
    epochs=30
    history = model.fit(
      train_ds,
      validation_data=val_ds,
      epochs=epochs,
      callbacks=[cp_callback]
    )

    ## Save Model ##
    model.save('savedmodels/fishNet')


    ## VISUALIZE RESULTS

    acc = history.history['accuracy']
    val_acc = history.history['val_accuracy']

    loss = history.history['loss']
    val_loss = history.history['val_loss']

    epochs_range = range(epochs)

    plt.figure(figsize=(8, 8))
    plt.subplot(1, 2, 1)
    plt.plot(epochs_range, acc, label='Training Accuracy')
    plt.plot(epochs_range, val_acc, label='Validation Accuracy')
    plt.legend(loc='lower right')
    plt.title('Training and Validation Accuracy')

    plt.subplot(1, 2, 2)
    plt.plot(epochs_range, loss, label='Training Loss')
    plt.plot(epochs_range, val_loss, label='Validation Loss')
    plt.legend(loc='upper right')
    plt.title('Training and Validation Loss')
    plt.show()

    return([class_names,model])


class_names = ["fish","nonfish"]
#[class_names,model] = train_model(batch_size,img_height,img_width)
print("Model training complete")



print(class_names)

## LOAD MODEL ##
model = tf.keras.models.load_model('savedmodels/fishNet')
model.summary()


## SLIDING WINDOW CLASSIFIER ##
from tensorflow.keras.preprocessing import image

# change this as you see fit
image_path = '/home/kostasl/workspace/zebrafishtrack/tensorDNN/img_target/00219.png'
#image_path = "/home/kostasl/workspace/zebrafishtrack/tensorDNN/img_target/templ_HB40_LR_camA_003_Templ_15752.jpg"

# Convert image to np.array
imgScene = image.load_img(image_path,color_mode = "grayscale")


#plt.figure(figsize=(6, 6))
#plt.imshow(image)
#plt.show()

#most deep learning models expect a batch of images as input
image_array = image.img_to_array(imgScene)
img_batch = np.expand_dims(image_array,axis=0)

imgb_size = np.shape(img_batch)
print(imgb_size)

#img_preprocessed = tf.keras.preprocessing.preprocess_input(img_batch)

# Sliding window
y_len,x_len,_ = image_array.shape

scale_x = img_width
scale_y = img_height


print("image size (%sx%s)",(y_len,x_len))
probability_model = Sequential([model,
                                layers.Softmax()])

## Sliding window of fixed size
for y in range(0, round(y_len/scale_y)*scale_y-scale_y,round(scale_y/2) ) :
    for x in range(0, round(x_len/scale_x)*scale_x-scale_x,round(scale_x/2) ):
        y_idxS = round(y)
        y_idxT = round(y+scale_y)
        x_idxS = round(x)
        x_idxT = round(x+scale_x)
        print('(%s,%s) -> (%s:%s),(%s:%s)' % (x + 1, y + 1,x_idxS,x_idxT,y_idxS,y_idxT))
        cropped_image = img_batch[:,y_idxS:y_idxT,
                                      x_idxS:x_idxT,:]
        #imgByteArray = io.BytesIO()
        #cropped_image.save(imgByteArray, format='JPEG')
        #imgByteArray = imgByteArray.getvalue()
        cimg_size = np.shape(cropped_image)
        #print(cimg_size)

        # Classify
        #classifier(imgByteArray,label_path,retrained_path)
        predictions = probability_model.predict(cropped_image)
        print(predictions[0])
        pclass = np.argmax(predictions[0])
        print("~~~ Predicted Class: %s %s" % (pclass,class_names[pclass] ) )

        if (pclass == 0):
            image.save_img("/home/kostasl/workspace/zebrafishtrack/tensorDNN/test/fish/00219-" + str(x)+"x"+str(y)+".jpg",cropped_image[0])
        else:
            image.save_img("/home/kostasl/workspace/zebrafishtrack/tensorDNN/test/nonfish/00219-" + str(x) + "x" + str(y) + ".jpg", cropped_image[0])




## END
input("Press Enter to exit...")



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
