# FishNet DNN - Train two tensorflow models :
#           * A model the fish anterior over all directions - to be used to locate the centre of fish template within a blob region
##          * B model of upright (Normed) fish antrerior - which is used to detect direction by testing against rotated image regions for best match
#import matplotlib.pyplot
#import io
#import PIL

import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib

from tensorflow.keras import datasets, layers, models
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras import regularizers
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.applications.efficientnet import EfficientNetB1, preprocess_input

# Templates 28x38 images size form the dataset, however data is augmented such that all fish rotations can be detected
# for this reason the training image needs to be square, so are rotations do not clip the top of the image (the eyes) when rotated
# 90 degrees. Altenrativelly, Classification could be done ny rotated the detected Regions and classify based on Normed - vertical orientiation only
#

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


# To Batch convert from pgm to jpg use mogrify :
# mogrify -format jpg *.pgm
#
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print('FishNet TensorFlow Model Training')
    print(tf.version.VERSION)

data_dir = pathlib.Path("/home/kostasl/workspace/zebrafishtrack/tensorDNN/trainset_cleaned/")
<<<<<<< HEAD
valid_dir = pathlib.Path("/home/kostasl/workspace/zebrafishtrack/tensorDNN/trainset/")
=======
valid_dir = pathlib.Path("/home/kostasl/workspace/zebrafishtrack/tensorDNN/trainset_cleaned/")
>>>>>>> 3bcea158f3182a7b751f5e736233034771dbee8f

fish = list(data_dir.glob('./fish/*.jpg'))
# PIL.Image.open(str(fish[0]))


nonfish = list(data_dir.glob('./nonfish/*.jpg'))
# PIL.Image.open(str(nonfish[0]))

bResetModelTraining = True  ## Do Not Incremental Train / Reset And Start over

batch_size = 32
img_height = 38
img_width = 28
epochs = 150
num_classes = 4
strModelPath = 'savedmodels/fishNet_loc_cleaned'
## Had To run x3 times with a validation split 0.3 - 0.5 before I got good filtering of entire scene - as tested by testModel
def train_model(epochs, batch_size, img_height, img_width, randRot=0.0
                , model=None):
    train_ds = tf.keras.preprocessing.image_dataset_from_directory(
        str(data_dir),
        validation_split=0.3,
        subset="training",
        seed=123,
        image_size=(img_height, img_width),
        color_mode='grayscale',
        batch_size=batch_size)

    val_ds = tf.keras.preprocessing.image_dataset_from_directory(
        str(data_dir),
        validation_split=0.3,
        subset="validation",
        seed=123,
        image_size=(img_height, img_width),
        color_mode='grayscale',
        batch_size=batch_size)
    class_names = train_ds.class_names

    train_datagen = tf.keras.preprocessing.image.ImageDataGenerator(horizontal_flip=True,
                                                                    validation_split=0.2)

    train_generator = train_datagen.flow_from_directory(
        str(data_dir),
        target_size=(img_height, img_width),
        color_mode='grayscale',
        batch_size=batch_size,
        class_mode='binary',
        classes=class_names)

    test_datagen = tf.keras.preprocessing.image.ImageDataGenerator()
    validation_generator = test_datagen.flow_from_directory(
        str(data_dir),
        target_size=(img_height, img_width),
        color_mode='grayscale',
        batch_size=batch_size,
        class_mode='binary',
        classes=class_names)

    y_labels = np.concatenate([y for x, y in train_ds], axis=0)
    x_images = np.concatenate([x for x, y in train_ds], axis=0)
    print(class_names)

    ##train_generator

    plt.figure(figsize=(10, 10))
    for images, labels in train_ds.take(1):
        for i in range(9):
            ax = plt.subplot(3, 3, i + 1)
            plt.imshow(images[i].numpy().astype("uint8"))  ##
            plt.title(class_names[labels[i]])
            plt.axis("off")

    for image_batch, labels_batch in train_ds:
        print(image_batch.shape)
        print(labels_batch.shape)
        break

    ## Data Augmentation
    ##Random Translation ##height, width representing lower and upper bound for shifting vertically.
    # A negative value means shifting image up, while a positive value means shifting image down.
    data_augmentation = keras.Sequential(
        [
            layers.experimental.preprocessing.RandomFlip("horizontal",
                                                         input_shape=(img_height, img_width, 1)),
            layers.experimental.preprocessing.RandomRotation(randRot),
            ##layers.experimental.preprocessing.RandomTranslation((0, 0.15), 0, fill_mode="nearest"),
            layers.experimental.preprocessing.RandomZoom(0.1),

        ]
    )

    fig = plt.figure(figsize=(10, 10))
    fig.suptitle("Augmented fish samples")
    ## Show Augmented Images of Fish samples ##
    for images, labels in train_ds.take(1):
        idxF = np.where(labels == 0)  ## Get All Idx Of all Fish Train samples
        for i in range(9):
            augmented_images = data_augmentation(images)
            ax = plt.subplot(3, 3, i + 1)
            plt.imshow(augmented_images[idxF[0][0]].numpy().astype("uint8"), label='Augmented data')
            plt.title(class_names[labels[idxF[0][0]]])
            plt.axis("off")

    # plt.show()

    ## Config Dataset performance
    AUTOTUNE = tf.data.AUTOTUNE

    train_ds = train_ds.cache().shuffle(1000).prefetch(buffer_size=AUTOTUNE)
    val_ds = val_ds.cache().prefetch(buffer_size=AUTOTUNE)

    ##Standardize the data
    normalization_layer = tf.keras.layers.experimental.preprocessing.Rescaling(1. / 255)

    ## Create the model ##

    if model is None:
        model = Sequential([
            data_augmentation,
            layers.experimental.preprocessing.Rescaling(1. / 255, input_shape=(img_height, img_width, 1)),
            layers.Conv2D(8, (3, 4), padding='same', activation='relu'),  # , activation='relu'
            #layers.Dropout(0.1),
            layers.MaxPooling2D(),
            layers.Conv2D(16, (3, 4), padding='same', activation='relu'),  # , activation='relu'
            #layers.Dropout(0.1),
            layers.MaxPooling2D(),
            layers.Conv2D(32, (3, 4), padding='same', activation='relu'),  # , activation='relu'
            #layers.Dropout(0.1),
            layers.MaxPooling2D(),
            layers.Conv2D(64, (3, 4), padding='same', activation='relu'),
            #layers.MaxPooling2D(),
            layers.Conv2D(128, (3, 4), padding='same', activation='relu'),
            layers.MaxPooling2D(),
            layers.Flatten(),
            layers.Dense(1250, activation='relu', kernel_regularizer=regularizers.l2(0.001)),
            layers.Dense(720, activation='relu', kernel_regularizer=regularizers.l2(0.001)),
            layers.Dropout(0.3),
            layers.Dense(200, activation='relu', kernel_regularizer=regularizers.l2(0.001)),
            layers.Dense(20, activation='relu', kernel_regularizer=regularizers.l2(0.001)),
            layers.Dropout(0.3),
            layers.Dense(num_classes, activation='softmax',name="output")  ##
        ])
        ##COMPILE MODEL
        model.compile(optimizer='adam',
                      # optimizer=RMSprop(lr=0.001),
                      #loss=tf.keras.losses.binary_crossentropy, ##For Binary Classifier
                      loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=False), ##Categorial classifier
                      metrics=['accuracy'])

    model.summary()

    checkpoint_path = "training_1/cp.ckpt"
    checkpoint_dir = os.path.dirname(checkpoint_path)

    # Create a callback that saves the model's weights
    cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath=checkpoint_path,
                                                     save_weights_only=True,
                                                     verbose=1)

    es_callback = EarlyStopping(monitor='val_loss',
                                patience=15)  ## Prevent Overfitting - Stop when Validation Loss drops for more than 5 epochs

    ## TRAIN MODEL
    history = model.fit(train_ds,
                        validation_data=val_ds,
                        epochs=epochs,
                        verbose=1,
                        callbacks=[cp_callback, es_callback]
                        )

    print(f'input_layer_name={model.input.name}')
    output_layer_name = model.output.name.split(':')[0]
    print(f'output_layer_name={output_layer_name}')

    # my_freeze_graph([output_layer_name], destination='savedmodels/frozen/', name="frozen_model.pb")

    ## VISUALIZE RESULTS
    acc = history.history['accuracy']
    val_acc = history.history['val_accuracy']

    loss = history.history['loss']
    val_loss = history.history['val_loss']

    epochs_range = range(np.shape(acc)[0])

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

    ## EVALUATE On DATASET CURVE ##
    #eval_results = model.evaluate(validation_generator)
    #print("EVALUATION: Test loss, Test acc:", eval_results)

    #validation_generator.reset()
    #validation_generator.shuffle()

    preds = model.predict(validation_generator, verbose=1)
    plt.subplot(1, 2, 1)
    plt.hist(preds, bins=2,label="Predictions Split")
    plt.subplot(1, 2, 2)
    plt.hist(validation_generator.classes, bins=2,label="Validation Data Split")

    #np.concatenate([y for x, y in train_ds], axis=0)
    #pEval = np.concatenate ([c for c in validation_generator.classes ] )
    #err = np.substract(pEval,preds)

    #plt.plot( err)
    #plt.show()

    return ([class_names, model])


def testModel(strTImg):
    ## Test SLIDING WINDOW CLASSIFIER ##
    from tensorflow.keras.preprocessing import image

    probability_model_loc = tf.keras.models.load_model(strModelPath)

    # change this as you see fit
    # image_path = pathlib.Path("/home/kostasl/workspace/hello_tf_c_api/build/nonfish_sample.jpg")
    image_path = pathlib.Path(strTImg)
    img_basename = os.path.basename(image_path)

    if (image_path.is_file()):
        # Convert image to np.array
        imgScene = image.load_img(image_path, color_mode="grayscale")
    else:
        print("ERROR image not found")
        raise ("Image not found")

    # plt.figure(figsize=(6, 6))
    # plt.imshow(image)
    # plt.show()

    # most deep learning models expect a batch of images as input
    image_array = image.img_to_array(imgScene)
    img_batch = np.expand_dims(image_array, axis=0)

    imgb_size = np.shape(img_batch)
    print(imgb_size)

    # img_preprocessed = tf.keras.preprocessing.preprocess_input(img_batch)

    # Sliding window
    y_len, x_len, _ = image_array.shape

    scale_x = img_width
    scale_y = img_height
    print("image size (%sx%s)", (y_len, x_len))

    ## Sliding window of fixed size
    bScanStop = False
    for y in range(0, round(y_len / scale_y) * scale_y, round(scale_y / 2)):
        if (bScanStop and ((y + scale_y) > y_len)):  ## Check If Boundary reached
            print("Image scan done")
            break
        else:
            bScanStop = False

        for x in range(0, round(x_len / scale_x) * scale_x, round(scale_x / 2)):
            if ((x + scale_x) > x_len):
                bScanStop = True
                break  ##Break Horizontal Scan Loop - reached boundary

            y_idxS = round(y)
            y_idxT = round(y + scale_y)
            x_idxS = round(x)
            x_idxT = round(x + scale_x)
            print('(%s,%s) -> (%s:%s),(%s:%s)' % (x + 1, y + 1, x_idxS, x_idxT, y_idxS, y_idxT))
            cropped_image = img_batch[:, y_idxS:y_idxT,
                            x_idxS:x_idxT, :]
            # imgByteArray = io.BytesIO()
            # cropped_image.save(imgByteArray, format='JPEG')
            # imgByteArray = imgByteArray.getvalue()
            cimg_size = np.shape(cropped_image)
            # print(cimg_size)

            # Classify
            # classifier(imgByteArray,label_path,retrained_path)
            predictions = probability_model_loc.predict(cropped_image)
            print(predictions[0])
            print(class_names)
            if (np.shape(predictions[0])[0] > 1):
                #pclass = np.argmax(predictions[0])
                if (predictions[0][0]-predictions[0][1] > 0.5):
                    pclass = 0
                else:
                    pclass = 1

                print("~~~ Predicted Class: %s %s" % (predictions[0], class_names[pclass]))
                image.save_img(
                "/home/kostasl/workspace/zebrafishtrack/tensorDNN/test/" + str(class_names[pclass]) + "/" + str(
                    img_basename) + "-" + str(x) + "x" + str(y) + ".jpg", cropped_image[0])
            else:
                if predictions[0] < 0.5:
                    print("~~~ Predicted Class: %s %s" % (predictions[0], class_names[0]))
                else:
                    print("~~~ Predicted Class: %s %s" % (predictions[0], class_names[1]))



    #            image.save_img("/home/kostasl/workspace/zebrafishtrack/tensorDNN/test/nonfish/00266-" + str(x) + "x" + str(y) + ".jpg", cropped_image[0])


class_names = ["fish", "nonfish"]


# # # #
model_dir_invar = None
model_directional = None
## LOAD MODEL ##
if (not bResetModelTraining):
    model_dir_invar = tf.keras.models.load_model(strModelPath)

## Train to identify fish in region regardless of orientation
[class_names, model_dir_invar] = train_model(epochs, batch_size, img_height, img_width, 1.0, model_dir_invar)
## Save Model ##
model_dir_invar.save(strModelPath)
# print("Model training complete")
print(class_names)

## GET OPERATION/LAYER NAMES INPUT -OUTPUT
# -Another way is by using saved_model_cli show --dir {mobilenet_save_path} --tag_set serve
if (not model_dir_invar is  None):
    model_dir_invar.summary()
    print(f'input_layer_name={model_dir_invar.input.name}')
    output_layer_name = model_dir_invar.output.name.split(':')[0]
    print(f'output_layer_name={output_layer_name},{model_dir_invar.output.name}')


print("~~~~~~~~~ Test Prediction on Non-fish and fish samples ~~~~")
testModel("/home/kostasl/workspace/zebrafishtrack/tensorDNN/valid/nonfish/00240-28x684.jpg")
testModel("/home/kostasl/workspace/zebrafishtrack/tensorDNN/valid/fish/templ_HB40_LR_camB_Templ_42695.jpg")
#testModel("/home/kostasl/workspace/zebrafishtrack/tensorDNN/valid/img_target/34450.png")
#testModel("/home/kostasl/workspace/zebrafishtrack/tensorDNN/valid/img_target/22271.png")
#testModel("/home/kostasl/workspace/zebrafishtrack/tensorDNN/valid/img_target/46047.png")


# ## MAKE PROB PREDICTION Version Add Softmax Layer And Save as fishNet_prob - THis version is used by the tracker
# print("Saving fishNet_prob Model With Probabilistic SOFTMAX output layer")
# probability_model_loc = Sequential([model_dir_invar,
#                                 layers.Softmax()])
# probability_model_loc.save('savedmodels/fishNet_loc_prob')
#
## ## DIRECTIONAL - UPRIGHT MODEL
## Save Model ##
#if (not bResetModelTraining):
# model_directional = tf.keras.models.load_model('savedmodels/fishNet_dir')
#
# ## Train to identify fish in region at a fixed (Vecrtical Orientation) Model used to Detect direction of fish
# [class_names,model_directional] = train_model(epochs,batch_size,img_height,img_width,0.0,model_directional)
# model_directional.save('savedmodels/fishNet_dir')
#
# ## MAKE PROB PREDICTION Version Add Softmax Layer And Save as fishNet_prob - THis version is used by the tracker
# # print("Saving fishNet_prob Directional Model With Probabilistic SOFTMAX output layer")
# # probability_model = Sequential([model_directional,
# #                                 layers.Softmax()])
# # probability_model.save('savedmodels/fishNet_dir_prob')
# if (not model_directional is  None):
#  print("Model training complete")
#  print(class_names)
#  model_directional.summary()
#  print(f'input_layer_name={model_directional.input.name}')
#  output_layer_name = model_directional.output.name.split(':')[0]
#  print(f'output_layer_name={output_layer_name},{model_directional.output.name}')
# #

## END
input("Press Enter to exit...")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
