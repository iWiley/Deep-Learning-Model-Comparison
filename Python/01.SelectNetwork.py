# This script is mainly used to filter the best model.
# The parameters related to the model can be viewed in real-time in TensorBoard.
import os
import pathlib

import pytorch_lightning as pl
import torch
import torch.nn as nn
from pytorch_lightning.callbacks import TQDMProgressBar
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from pytorch_lightning.loggers import TensorBoardLogger
from torch.utils.data import DataLoader, WeightedRandomSampler
from torchvision import transforms
from torchvision.datasets import ImageFolder

from Python.script.Model import Model

dataDir = './Python/data/Train'
if len(os.listdir(dataDir)) == 0:
    raise Exception(
        "Please ensure that the data directory 'Python/data/Train' is not empty. You may need to download training data from the Github release page of this project.")

print("The default parameters may not be suitable for your computer. If you cannot execute them successfully, please modify these parameters.")
batch_size = 150
batch_size_Inception_V3 = 20
num_workers = 12
dorpOutValue = .5
opt_func = torch.optim.AdamW
lr = .0001
# Ensure that the results are reproducible.
pl.seed_everything(666666)
# default image size except Inception_V3
tran_150 = transform = transforms.Compose([
    transforms.Resize(
        (150, 150), interpolation=transforms.InterpolationMode.BICUBIC),
    transforms.ToTensor()
])
# The minimum size required for the Inception_V3 is 299
tran_299 = transform = transforms.Compose([
    transforms.Resize(
        (299, 299), interpolation=transforms.InterpolationMode.BICUBIC),
    transforms.ToTensor()
])
# Split Dataset
full_dataset = ImageFolder(dataDir)
train_size = int(0.8 * len(full_dataset))
test_size = len(full_dataset) - train_size
train_ds, test_ds = torch.utils.data.random_split(
    full_dataset, [train_size, test_size])
# Number of classes
root = pathlib.Path(dataDir)
classes = sorted([j.name.split('/')[-1] for j in root.iterdir()])
# Here is the number of each classes of the training set. Because it is time-consuming, the source code has been commented and replaced with the execution result.
# train_classes = [label for _, label in train_ds]
# class_count = Counter(train_classes)
# class_weights = torch.Tensor(
# print(class_weights)
class_weights = 1 / torch.Tensor([58840, 66217, 7039, 138])
sample_weights = [0] * len(train_ds)
print("Calculating weight...")
# It is also extremely slow because it needs to traverse every sample.
for idx, (image, label) in enumerate(train_ds):
    class_weight = class_weights[label]
    sample_weights[idx] = class_weight
sampler = WeightedRandomSampler(weights=sample_weights,
                                num_samples=len(train_ds), replacement=True)
train_dl = DataLoader(train_ds, batch_size,
                      num_workers=num_workers, pin_memory=True, sampler=sampler)
val_dl = DataLoader(
    test_ds, batch_size, num_workers=num_workers, pin_memory=True)
# All models to be tested.
models = [
    Model.Model_Inception_V3.Model_VGG19,
    Model.Model_ResNet_18,
    Model.Model_ResNet_50,
    Model.Model_ResNet_101,
    Model.Model_ResNet_152,
    Model.Model_ModelNet_V3,
    Model.Model_ShuffleNet_V2,
    Model.Model_Wide_ResNet50_2,
    Model.Model_Wide_ResNet101_2,
    Model.Model_Inception_V3,
]


# To facilitate observation of training progress.
class ProgressBarEx(TQDMProgressBar):
    def __init__(self, s: str):
        super().__init__(10, 0)
        self._s = s

    def init_train_tqdm(self):
        bar = super().init_train_tqdm()
        bar.set_description(f'Training:{self._s}')
        return bar

    def init_validation_tqdm(self):
        bar = super().init_validation_tqdm()
        bar.set_description(f'Validating:{self._s}')
        return bar


# Start training
print("training...")
for m in models:
    trainningModel = m(len(classes), dorpOutValue, opt_func, lr)
    full_dataset.transform = \
        tran_299 if trainningModel.ModelName() == 'InceptionV3' else tran_150
    if trainningModel.ModelName() == 'InceptionV3':
        train_dl = DataLoader(train_ds, batch_size_Inception_V3,
                              num_workers=num_workers, pin_memory=True, sampler=sampler)
        val_dl = DataLoader(
            test_ds, batch_size_Inception_V3, num_workers=num_workers, pin_memory=True)
    trainningModel.loss = nn.CrossEntropyLoss(class_weights)
    logger = TensorBoardLogger(
        'TensorBoardLogger', prefix=trainningModel.ModelName())
    trainer = pl.Trainer(
        devices=-1, accelerator="auto", logger=logger,
        callbacks=[
            EarlyStopping(monitor="valid_loss_epoch",
                          patience=5,
                          min_delta=.0005,
                          check_on_train_epoch_end=True),
            ProgressBarEx(trainningModel.ModelName())])
    trainer.fit(trainningModel, train_dataloaders=train_dl,
                val_dataloaders=val_dl)
    trainer.save_checkpoint(
        f'./Python/models/{trainningModel.ModelName()}.pth')
