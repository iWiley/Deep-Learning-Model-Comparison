# Estimate the immune infiltration of WSIs in TCGA cohort.

import csv
import os

import numpy as np
import pytorch_lightning as pl
import torch
from torchvision import transforms

from Python.script.Model import Model_ShuffleNet_V2
from Python.script.TCGAWSILoadHelper import TCGAWSILoadHelper

wsipath = 'Python/data/TCGA'
wsis = os.listdir(wsipath)
wsis = [wsipath + x for x in wsis if isinstance(x, str)]
if len(wsis) == 0:
    raise Exception(
        "Please ensure that the data directory 'Python/data/TCGA' is not empty. You may need to download training data from the Github release page of this project.")

print("The default parameters may not be suitable for your computer. If you cannot execute them successfully, please modify these parameters.")
batch_size = 50
num_workers = 12
pl.seed_everything(666666)
trainer = pl.Trainer(devices=-1, accelerator="auto")
torch.set_float32_matmul_precision("medium")
ws = TCGAWSILoadHelper(wsis, batch_size, num_workers, None)
transform = transforms.Compose([
    transforms.Resize(
        (150, 150), interpolation=transforms.InterpolationMode.BICUBIC),
    transforms.ToTensor(),
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
ws.batch_size = batch_size
ws.trans = transform
# Check save location.
if os.path.exists('Python/score/TCGA.csv'):
    os.remove('Python/score/TCGA.csv')
file = open('Python/score/TCGA.csv', 'w')
writer = csv.writer(file)
writer.writerow(['File', 'invalid', 'tumor', 'TIL', 'TLS'])
# Start evaluation.
for dl in ws:
    print(ws.count, len(ws.wsis))
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
    writer.writerow([ws.WSIFile, v0, v1, v2, v3])
    file.flush()
