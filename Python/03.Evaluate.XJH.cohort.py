# Estimate the immune infiltration of WSIs in XJH cohort.

import csv
import os

import numpy as np
import pytorch_lightning as pl
import torch
from torch.utils.data import DataLoader
from torchvision import transforms

from Python.script.LocalDS import LocalDS
from Python.script.Model import Model_ShuffleNet_V2

wsipath = 'Python/data/XJH'
wsis = os.listdir(wsipath)
wsis = [wsipath + '/' + x for x in wsis if isinstance(x, str)]
wsis = [x for x in wsis if os.path.isdir(x)]
if len(wsis) == 0:
    raise Exception(
        "Please ensure that the data directory 'Python/data/XJH' is not empty. You may need to download training data from the Github release page of this project.")

pl.seed_everything(666666)
batch_size = 50
batch_size_Inception_V3 = 20
num_workers = 12

trainer = pl.Trainer(devices=-1, accelerator="auto")
torch.set_float32_matmul_precision("medium")
transform = transforms.Compose([
    transforms.Resize(
        (150, 150), interpolation=transforms.InterpolationMode.BICUBIC),
    transforms.ToTensor()
])
# Initialize default parameters and load the model.
model = Model_ShuffleNet_V2(4, .5, torch.optim.AdamW, .0001)
checkpoint = torch.load('Python/models/ShuffleNet_V2.pth')
# Because loss.weight is not a fixed element in the model,
# it cannot be loaded and needs to be removed.
# loss.weight only works in training.
# Removing it here does not affect the disease evaluation of WSIs
m = {k: v for k, v in checkpoint['state_dict'].items() if k != 'loss.weight'}
model.load_state_dict(m)
# Check save location.
if os.path.exists('Python/score/XJH.csv'):
    os.remove('Python/score/XJH.csv')
file = open('Python/score/XJH.csv', 'w')
writer = csv.writer(file)
writer.writerow(['File', 'invalid', 'tumor', 'TIL', 'TLS'])
# Start evaluation.
for wsfile in wsis:
    full_dataset = LocalDS(wsfile, transform=transform)
    dl = DataLoader(full_dataset, batch_size,
                    num_workers=num_workers, pin_memory=True)
    result = trainer.predict(model, dataloaders=dl)
    re = []
    for r in result:
        re.extend(r.numpy())
    re = np.argmax(re, axis=1)
    v0 = 0
    v1 = 0
    v2 = 0
    v3 = 0
    for v in re:
        if v == 0:
            v0 += 1
        elif v == 1:
            v1 += 1
        elif v == 2:
            v2 += 1
        else:
            v3 += 1
    writer.writerow([wsfile, v0, v1, v2, v3])
    file.flush()
