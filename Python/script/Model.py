# Various model networks and specific training steps.
from abc import abstractmethod

import numpy as np
import pytorch_lightning as pl
import torch.nn as nn
import torchmetrics
from Draw import DrawPRR, DrawROC
from torchmetrics import F1Score, PrecisionRecallCurve
from torchmetrics.classification import MulticlassConfusionMatrix
from torchvision import models


class Model(pl.LightningModule):
    loss = nn.CrossEntropyLoss()

    def __init__(self, num_classes, opt_func, lr):
        super().__init__()
        self.opt_func = opt_func
        self.lr = lr
        self.loss = nn.CrossEntropyLoss()
        self.train_acc = torchmetrics.Accuracy(
            task='multiclass', num_classes=num_classes)
        self.valid_acc = torchmetrics.Accuracy(
            task='multiclass', num_classes=num_classes)
        self.roc = torchmetrics.ROC(
            task='multiclass', num_classes=num_classes)
        self.auc = torchmetrics.AUROC(
            task='multiclass', num_classes=num_classes, average=None)
        self.mcm = MulticlassConfusionMatrix(
            num_classes=num_classes)
        self.mlprc = PrecisionRecallCurve(
            task='multiclass', num_classes=num_classes)
        self.f1 = F1Score(task="multiclass", num_classes=num_classes)

    @abstractmethod
    def ModelName(self) -> str:
        return ''

    def forward(self, xb):
        return self.network(xb)

    def training_step(self, batch, batch_idx):
        images, labels = batch
        out = self(images)
        if (self.ModelName() == 'InceptionV3'):
            out = out[0]
        loss = self.loss(out, labels)
        self.train_acc(out, labels)
        self.log('train_loss', loss, on_step=True,
                 on_epoch=True, prog_bar=True, logger=True)
        self.log('train_acc', self.train_acc, on_step=True,
                 on_epoch=True, prog_bar=True, logger=True)
        return loss

    def training_epoch_end(self, outputs):
        self.log('train_acc_epoch', self.train_acc.compute(),
                 on_epoch=True, logger=True)
        self.train_acc.reset()

    def validation_step(self, batch, batch_idx):
        images, labels = batch
        out = self(images)
        loss = self.loss(out, labels)
        self.valid_acc(out, labels)
        self.roc(out, labels)
        self.auc(out, labels)
        self.mcm(out, labels)
        self.mlprc(out, labels)
        self.f1(out, labels)
        self.log('valid_acc', self.valid_acc, on_step=True,
                 on_epoch=True, prog_bar=True, logger=True)
        self.log("valid_loss", loss, on_step=True,
                 on_epoch=True, prog_bar=True, logger=True)
        self.log("valid_f1", self.f1, on_step=True,
                 on_epoch=True, prog_bar=True, logger=True)

    def validation_epoch_end(self, outputs):
        self.valid_acc.reset()
        self.f1.reset()
        fpr, tpr, _ = self.roc.compute()
        auc = self.auc.compute()
        f = DrawROC(fpr, tpr, auc)
        f.savefig(
            f'svgs/{self.ModelName()}-{self.current_epoch}.ROC.svg',
            format="svg",
            dpi=300)

        self.roc.reset()
        self.auc.reset()

        mcm = self.mcm.compute()
        Q = np.array(mcm.detach().numpy())
        np.savetxt(f'svgs/{self.ModelName()}-{self.current_epoch}.MCM.txt', Q)

        precision, recall, _ = self.mlprc.compute()
        f = DrawPRR(precision, recall)
        f.savefig(
            f'svgs/{self.ModelName()}-{self.current_epoch}.PRR.svg',
            format="svg",
            dpi=300)
        self.mcm.reset()
        self.mlprc.reset()

    def get_progress_bar_dict(self):
        tqdm_dict = super().get_progress_bar_dict()
        tqdm_dict.pop("loss", None)
        tqdm_dict.pop("v_num", None)
        return tqdm_dict

    def configure_optimizers(self):
        return self.opt_func(self.parameters(), lr=self.lr)


def _SetGrad(model: nn.Module):
    for p in model.parameters():
        p.requires_grad = False


class Model_VGG19(Model):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model_vgg = models.vgg19(weights=models.VGG19_Weights.DEFAULT)
        _SetGrad(model_vgg)
        model_vgg.classifier = nn.Sequential(
            nn.Linear(in_features=25088, out_features=2048),
            nn.ReLU(),
            nn.Linear(in_features=2048, out_features=512),
            nn.ReLU(),
            nn.Dropout(p=dropOutValue),
            nn.Linear(in_features=512, out_features=num_classes),
            nn.LogSoftmax(dim=1)
        )
        super().__init__(num_classes, opt_func, lr)
        self.network = model_vgg

    def ModelName(self) -> str:
        return "VGG19"


class Model_ResNet(Model):
    def __init__(self, model: nn.Module, num_classes, dropOutValue, opt_func, lr):
        _SetGrad(model)
        model.fc = nn.Sequential(
            nn.Linear(in_features=2048, out_features=1024),
            nn.ReLU(),
            nn.Linear(in_features=1024, out_features=512),
            nn.ReLU(),
            nn.Dropout(p=dropOutValue),
            nn.Linear(in_features=512, out_features=num_classes),
            nn.LogSoftmax(dim=1))
        super().__init__(num_classes, opt_func, lr)
        self.network = model


class Model_ResNet_18(Model):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.resnet18(weights=models.ResNet18_Weights.DEFAULT)
        print(model)
        _SetGrad(model)
        model.fc = nn.Sequential(
            nn.Linear(in_features=512, out_features=128),
            nn.ReLU(),
            nn.Dropout(p=dropOutValue),
            nn.Linear(in_features=128, out_features=num_classes),
            nn.LogSoftmax(dim=1))
        super().__init__(num_classes, opt_func, lr)
        self.network = model

    def ModelName(self) -> str:
        return "ResNet18"


class Model_ResNet_50(Model_ResNet):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.resnet50(weights=models.ResNet50_Weights.DEFAULT)
        super().__init__(model, num_classes, dropOutValue, opt_func, lr)

    def ModelName(self) -> str:
        return "ResNet50"


class Model_ResNet_101(Model_ResNet):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.resnet101(weights=models.ResNet101_Weights.DEFAULT)
        super().__init__(model, num_classes, dropOutValue, opt_func, lr)

    def ModelName(self) -> str:
        return "ResNet101"


class Model_ResNet_152(Model_ResNet):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.resnet152(weights=models.ResNet152_Weights.DEFAULT)
        super().__init__(model, num_classes, dropOutValue, opt_func, lr)

    def ModelName(self) -> str:
        return "ResNet152"


class Model_ModelNet_V3(Model):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.mobilenet_v3_large(
            weights=models.MobileNet_V3_Large_Weights.DEFAULT)
        _SetGrad(model)
        model.classifier = nn.Sequential(
            nn.Linear(in_features=960, out_features=512, bias=True),
            nn.Hardswish(),
            nn.Dropout(p=dropOutValue),
            nn.Linear(in_features=512, out_features=num_classes),
            nn.LogSoftmax(dim=1)
        )
        super().__init__(num_classes, opt_func, lr)
        self.network = model

    def ModelName(self) -> str:
        return "ModelNetV3"


class Model_Inception_V3(Model):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.inception_v3(
            weights=models.Inception_V3_Weights.DEFAULT)
        _SetGrad(model)
        model.fc = nn.Sequential(
            nn.Linear(in_features=2048, out_features=1024),
            nn.ReLU(),
            nn.Linear(in_features=1024, out_features=512),
            nn.ReLU(),
            nn.Dropout(p=dropOutValue),
            nn.Linear(in_features=512, out_features=num_classes),
            nn.LogSoftmax(dim=1))
        super().__init__(num_classes, opt_func, lr)
        self.network = model

    def ModelName(self) -> str:
        return "InceptionV3"


class Model_Wide_ResNet50_2(Model):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.wide_resnet50_2()
        _SetGrad(model)
        model.fc = nn.Sequential(
            nn.Linear(in_features=2048, out_features=1024),
            nn.ReLU(),
            nn.Linear(in_features=1024, out_features=512),
            nn.ReLU(),
            nn.Dropout(p=dropOutValue),
            nn.Linear(in_features=512, out_features=num_classes),
            nn.LogSoftmax(dim=1))
        super().__init__(num_classes, opt_func, lr)
        self.network = model

    def ModelName(self) -> str:
        return "Wide_ResNet50_2"


class Model_Wide_ResNet101_2(Model):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.wide_resnet101_2()
        _SetGrad(model)
        model.fc = nn.Sequential(
            nn.Linear(in_features=2048, out_features=1024),
            nn.ReLU(),
            nn.Linear(in_features=1024, out_features=512),
            nn.ReLU(),
            nn.Dropout(p=dropOutValue),
            nn.Linear(in_features=512, out_features=num_classes),
            nn.LogSoftmax(dim=1))
        super().__init__(num_classes, opt_func, lr)
        self.network = model

    def ModelName(self) -> str:
        return "Wide_ResNet101_2"


class Model_ShuffleNet_V2(Model):
    def __init__(self, num_classes, dropOutValue, opt_func, lr):
        model = models.shufflenet_v2_x1_0(
            weights=models.ShuffleNet_V2_X1_0_Weights.DEFAULT)
        print(model)
        _SetGrad(model)
        model.fc = nn.Sequential(
            nn.Linear(in_features=1024, out_features=50),
            nn.ReLU(),
            nn.Dropout(p=dropOutValue),
            nn.Linear(in_features=50, out_features=num_classes),
            nn.LogSoftmax(dim=1))
        super().__init__(num_classes, opt_func, lr)
        self.network = model

    def ModelName(self) -> str:
        return "ShuffleNet_V2"
