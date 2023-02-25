# To draw ROC curve and PR curve, 
# please uncomment relevant code in Model.py to generate pictures.
import matplotlib.pyplot as plt

def plot_AUC(ax, FPR, TPR, AUC, label):
	label = f'{str(label)}: {AUC:.3f}'
	ax.plot(FPR, TPR, label = label)
	return ax

def DrawROC(FPR, TPR, AUC):
    plt.style.use('ggplot')
    fig, ax = plt.subplots()

    for i in range(len(FPR)):
        ax = plot_AUC(ax, FPR[i].cpu().numpy(), TPR[i].cpu().numpy(), AUC[i].cpu().numpy(), label = str(i))

    ax.plot((0, 1), (0, 1), ':', color = 'grey')
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_aspect('equal')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend()
    return fig

def plot_PRR(ax, precision, recall):
	ax.plot(precision, recall)
	return ax

def DrawPRR(precision, recall):
    plt.style.use('ggplot')
    fig, ax = plt.subplots()

    for i in range(len(precision)):
        ax = plot_AUC(ax, precision[i].cpu().numpy(), recall[i].cpu().numpy())

    ax.plot((0, 1), (0, 1), ':', color = 'grey')
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_aspect('equal')
    ax.set_xlabel('Precision')
    ax.set_ylabel('Recall')
    ax.legend()
    return fig