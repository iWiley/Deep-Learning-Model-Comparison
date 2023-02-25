# Used to load data locally
import os

from PIL import Image
from torch.utils.data import Dataset


class LocalDS(Dataset):
    def __init__(self, pth, transform):
        super().__init__()
        wsis = os.listdir(pth)
        wsis = [pth + '/' + x for x in wsis if isinstance(x, str)]
        self.imgs = [x for x in wsis if os.path.isfile(x)]
        self.tramsfrom = transform

    def __getitem__(self, index):
        img = Image.open(self.imgs[index])
        img = img.convert("RGB")
        img = self.tramsfrom(img)
        return img

    def __len__(self):
        return len(self.imgs)
