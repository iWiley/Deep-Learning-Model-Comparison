# Used to automatically return TCGAWSIDataset to simplify operation
from threading import Thread

import TCGAWSIDataset
from torch.utils.data import DataLoader


class TCGAWSILoadHelper:
    def __init__(self, wsis, batch_size, num_workers, trans):
        self.wsis = wsis
        self.count = 0
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.trans = trans
        self.t = None
        self.__load_new__()

    def __iter__(self):
        return self

    def __next__(self):
        if self.count < len(self.wsis):
            if (self.t != None):
                self.t.join()
            self.t = Thread(target=self.__load_new__())
            self.count += 1
            self.t.start()
            return self.current
        else:
            raise StopIteration

    def __load_new__(self):
        ds = TCGAWSIDataset(
            wsiFile=self.wsis[self.count],
            transform=self.trans)
        self.current = DataLoader(
            ds, self.batch_size, num_workers=self.num_workers, pin_memory=True)
        self.t = None

    @property
    def WSIFile(self):
        return self.wsis[self.count-1]
