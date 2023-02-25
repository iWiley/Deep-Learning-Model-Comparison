# Used to load data from a WSI file.
import openslide
from torch.utils.data import Dataset
import math
import torchvision.transforms as transforms

toTensor = transforms.ToTensor()


class TCGAWSIDataset(Dataset):
    _tileSize = (512, 512)
    _layer = 1
    _zoom: int
    slide: openslide.OpenSlide
    xcount: int
    ycount: int

    def __init__(self, wsiFile: str, transform):
        super().__init__()
        self.slide = openslide.OpenSlide(wsiFile)
        xtotal = self.slide.level_dimensions[1][0]
        ytotal = self.slide.level_dimensions[1][1]
        self._zoom = self.slide.level_downsamples[1] / \
            self.slide.level_downsamples[0]
        self.xcount = int(math.floor(float(xtotal) / (self._tileSize[0])))
        self.ycount = int(math.floor(float(ytotal) / (self._tileSize[1])))
        self.transform = transform

    def __len__(self):
        return self.xcount * self.ycount

    def __getitem__(self, idx):
        x = int((idx % self.xcount) * self._tileSize[0] * self._zoom)
        y = int(math.floor(idx / self.xcount) * self._tileSize[1] * self._zoom)
        item = self.slide.read_region((x, y), self._layer, self._tileSize)
        item = item.convert("RGB")
        if self.transform is not None:
            item = self.transform(item)
        return item
