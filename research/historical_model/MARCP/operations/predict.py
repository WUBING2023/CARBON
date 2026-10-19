"""Prediction from sequence: one-head and two-head (treat/control), fragment or midBase mode."""

import torch
import numpy as np

from MARCP.utils.functions import onehot_encode, sequence2int


def predict2HeadW(seqTranscribed, window_size, signalList, model, mode, codeModel="OneHot", device=None):
    remainder = (len(seqTranscribed) % window_size)
    print(f"## remainder Length is: {remainder}")

    if mode == "Fragment":
        num_sublists = len(seqTranscribed) // window_size
        sublists = [seqTranscribed[i * window_size: (i + 1) * window_size] for i in range(num_sublists)]
        remainder = len(seqTranscribed) % window_size
        if remainder:
            sublists.append(seqTranscribed[-window_size:])

    elif mode == "midBase":
        sublists = [seqTranscribed[i: i + window_size] for i in np.arange(0, len(seqTranscribed) - window_size + 1, step=1)]

    if codeModel == "OneHot":
        vocab = list("AUGC")
        seqsInt = onehot_encode(sublists, vocab, 4)
        subseqOH = torch.from_numpy(np.asarray(seqsInt)).to(device)
    elif codeModel == "Embedding":
        seqsInt = [sequence2int(seq) for seq in sublists]
    subseqOH = torch.from_numpy(np.asarray(seqsInt)).to(device)

    model.eval()
    with torch.no_grad():
        temp = []
        if (codeModel == "Embedding"):
            subseqOH = subseqOH.int()
        elif (codeModel == "OneHot"):
            subseqOH = subseqOH.float()

        if mode == "Fragment":
            if remainder:
                for subseqOHTurn in subseqOH[:-1].split(128, 0):
                    outputs = model(subseqOHTurn)
                    w = outputs[2]
                    w = w.unsqueeze(axis=1)
                    treat_predict = (outputs[0] * w) + (1 - w) * outputs[1]
                    treat_predict = treat_predict.detach().cpu().numpy().tolist()
                    temp.extend(treat_predict)
                temp = np.array(temp).reshape(-1).tolist()
                print("first temp len", len(temp))
                outputs_end = model(subseqOH[-1:])
                w_end = outputs_end[2]
                w_end = w_end.unsqueeze(axis=1)
                treat_predict_end = (outputs_end[0] * w_end) + (1 - w_end) * outputs_end[1]
                treat_predict_end = treat_predict_end.detach().cpu().numpy().tolist()
                print(remainder, len(temp), len(np.array(treat_predict_end).reshape(-1).tolist()[-remainder:]), len(np.array(treat_predict_end).reshape(-1).tolist()))
                temp = temp + np.array(treat_predict_end).reshape(-1).tolist()[-remainder:]
                print("remainder: ", remainder, len(treat_predict_end[-remainder:]), len(temp))
            else:
                for subseqOHTurn in subseqOH.split(128, 0):
                    outputs = model(subseqOHTurn)
                    w = outputs[2]
                    w = w.unsqueeze(axis=1)
                    treat_predict = (outputs[0] * w) + (1 - w) * outputs[1]
                    treat_predict = treat_predict.detach().cpu().numpy().tolist()
                    temp.extend(treat_predict)
                temp = np.array(temp).reshape(-1).tolist()

        elif mode == "midBase":
            for subseqOHTurn in subseqOH.split(1000, 0):
                subseqPredicted = model(subseqOHTurn)
                outputs = subseqPredicted
                w = outputs[2]
                w = w.unsqueeze(axis=1)
                treat_predict = (outputs[0] * w) + (1 - w) * outputs[1]
                treat_predict = treat_predict[:, window_size // 2]
                treat_predict = treat_predict.detach().cpu().numpy().tolist()
                temp.extend(treat_predict)
            remainder_index = 0
            outputs_end = model(subseqOH[remainder_index].unsqueeze(0))
            w_end = outputs_end[2]
            w_end = w_end.unsqueeze(axis=1)
            treat_predict_end = (outputs_end[0] * w_end) + (1 - w_end) * outputs_end[1]
            treat_predict_end = treat_predict_end.detach().cpu().numpy()
            temp = treat_predict_end[:, :window_size // 2].reshape(-1).tolist() + temp
            remainder_index = -1
            outputs_end = model(subseqOH[remainder_index].unsqueeze(0))
            w_end = outputs_end[2]
            w_end = w_end.unsqueeze(axis=1)
            treat_predict_end = (outputs_end[0] * w_end) + (1 - w_end) * outputs_end[1]
            treat_predict_end = treat_predict_end.detach().cpu().numpy()
            temp = treat_predict_end[:, :window_size // 2].reshape(-1).tolist() + temp

        transcriptPd = np.array(temp)
        print("Signal lenth of prediction and transcript length !", len(transcriptPd), len(signalList))
        if (len(transcriptPd) != len(signalList)):
            print("Error: signal lenth not equal transcript length !", len(transcriptPd), len(signalList))
            return False
    return transcriptPd


def predict1Head(seqTranscribed, window_size, signalList, model, mode, codeModel="OneHot", device=None, batch=10):
    remainder = (len(seqTranscribed) % window_size)
    print(f"## remainder Length is: {remainder}")

    if mode == "Fragment":
        num_sublists = len(seqTranscribed) // window_size
        sublists = [seqTranscribed[i * window_size: (i + 1) * window_size] for i in range(num_sublists)]
        remainder = len(seqTranscribed) % window_size
        if remainder:
            sublists.append(seqTranscribed[-window_size:])

    elif mode == "midBase":
        sublists = [seqTranscribed[i: i + window_size] for i in np.arange(0, len(seqTranscribed) - window_size + 1, step=1)]

    if codeModel == "OneHot":
        vocab = list("AUGC")
        seqsInt = onehot_encode(sublists, vocab, 4)
        subseqOH = torch.from_numpy(np.asarray(seqsInt)).to(device)
    elif codeModel == "Embedding":
        seqsInt = [sequence2int(seq) for seq in sublists]
    subseqOH = torch.from_numpy(np.asarray(seqsInt)).to(device)

    model.eval()
    with torch.no_grad():
        temp = []
        if (codeModel == "Embedding"):
            subseqOH = subseqOH.int()
        elif (codeModel == "OneHot"):
            subseqOH = subseqOH.float()

        if mode == "Fragment":
            if remainder:
                for subseqOHTurn in subseqOH[:-1].split(batch, 0):
                    treat_predict = model(subseqOHTurn)
                    treat_predict = treat_predict.detach().cpu().numpy().tolist()
                    temp.extend(treat_predict)
                temp = np.array(temp).reshape(-1).tolist()
                treat_predict_end = model(subseqOH[-1:])
                treat_predict_end = treat_predict_end.detach().cpu().numpy().tolist()
                print(remainder, len(temp), len(np.array(treat_predict_end).reshape(-1).tolist()[-remainder:]), len(np.array(treat_predict_end).reshape(-1).tolist()))
                temp = temp + np.array(treat_predict_end).reshape(-1).tolist()[-remainder:]
                print("remainder: ", remainder, len(treat_predict_end[-remainder:]), len(temp))
            else:
                for subseqOHTurn in subseqOH.split(1000, 0):
                    treat_predict = model(subseqOHTurn)
                    treat_predict = treat_predict.detach().cpu().numpy().tolist()
                    temp.extend(treat_predict)
                temp = np.array(temp).reshape(-1).tolist()

        elif mode == "midBase":
            for subseqOHTurn in subseqOH.split(batch, 0):
                treat_predict = model(subseqOHTurn)
                treat_predict = treat_predict[:, window_size // 2]
                treat_predict = treat_predict.detach().cpu().numpy().tolist()
                temp.extend(treat_predict)
            remainder_index = 0
            treat_predict_end = model(subseqOH[remainder_index].unsqueeze(0))
            treat_predict_end = treat_predict_end.detach().cpu().numpy()
            temp = treat_predict_end[:, :window_size // 2].reshape(-1).tolist() + temp
            remainder_index = -1
            treat_predict_end = model(subseqOH[remainder_index].unsqueeze(0))
            treat_predict_end = treat_predict_end.detach().cpu().numpy()
            temp = treat_predict_end[:, window_size // 2:].reshape(-1).tolist() + temp

        transcriptPd = np.array(temp)
        print("Signal lenth of prediction and transcript length !", len(transcriptPd), len(signalList))
        if (len(transcriptPd) != len(signalList)):
            print("Error: signal lenth not equal transcript length !", len(transcriptPd), len(signalList))
            return False
    return transcriptPd


def predict2HeadWFromSeq(seqTranscribed, window_size, model, mode, codeModel="OneHot", device=None):
    remainder = (len(seqTranscribed) % window_size)
    print(f"## remainder Length is: {remainder}")

    if mode == "Fragment":
        num_sublists = len(seqTranscribed) // window_size
        sublists = [seqTranscribed[i * window_size: (i + 1) * window_size] for i in range(num_sublists)]
        remainder = len(seqTranscribed) % window_size
        if remainder:
            sublists.append(seqTranscribed[-window_size:])

    elif mode == "midBase":
        sublists = [seqTranscribed[i: i + window_size] for i in np.arange(0, len(seqTranscribed) - window_size + 1, step=1)]

    if codeModel == "OneHot":
        vocab = list("AUGC")
        seqsInt = onehot_encode(sublists, vocab, 4)
        subseqOH = torch.from_numpy(np.asarray(seqsInt)).to(device)
    elif codeModel == "Embedding":
        seqsInt = [sequence2int(seq) for seq in sublists]
    subseqOH = torch.from_numpy(np.asarray(seqsInt)).to(device)

    model.eval()
    with torch.no_grad():
        temp = []
        if (codeModel == "Embedding"):
            subseqOH = subseqOH.int()
        elif (codeModel == "OneHot"):
            subseqOH = subseqOH.float()

        if mode == "Fragment":
            if remainder:
                for subseqOHTurn in subseqOH[:-1].split(1000, 0):
                    outputs = model(subseqOHTurn)
                    w = outputs[2]
                    w = w.unsqueeze(axis=1)
                    treat_predict = (outputs[0] * w) + (1 - w) * outputs[1]
                    treat_predict = treat_predict.detach().cpu().numpy().tolist()
                    temp.extend(treat_predict)
                temp = np.array(temp).reshape(-1).tolist()
                print("first temp len", len(temp))
                outputs_end = model(subseqOH[-1:])
                w_end = outputs_end[2]
                w_end = w_end.unsqueeze(axis=1)
                treat_predict_end = (outputs_end[0] * w_end) + (1 - w_end) * outputs_end[1]
                treat_predict_end = treat_predict_end.detach().cpu().numpy().tolist()
                print(remainder, len(temp), len(np.array(treat_predict_end).reshape(-1).tolist()[-remainder:]), len(np.array(treat_predict_end).reshape(-1).tolist()))
                temp = temp + np.array(treat_predict_end).reshape(-1).tolist()[-remainder:]
                print("remainder: ", remainder, len(treat_predict_end[-remainder:]), len(temp))
            else:
                for subseqOHTurn in subseqOH.split(1000, 0):
                    outputs = model(subseqOHTurn)
                    w = outputs[2]
                    w = w.unsqueeze(axis=1)
                    treat_predict = (outputs[0] * w) + (1 - w) * outputs[1]
                    treat_predict = treat_predict.detach().cpu().numpy().tolist()
                    temp.extend(treat_predict)
                temp = np.array(temp).reshape(-1).tolist()

        elif mode == "midBase":
            for subseqOHTurn in subseqOH.split(1000, 0):
                subseqPredicted = model(subseqOHTurn)
                outputs = subseqPredicted
                w = outputs[2]
                w = w.unsqueeze(axis=1)
                treat_predict = (outputs[0] * w) + (1 - w) * outputs[1]
                treat_predict = treat_predict[:, window_size // 2]
                treat_predict = treat_predict.detach().cpu().numpy().tolist()
                temp.extend(treat_predict)
            remainder_index = 0
            outputs_end = model(subseqOH[remainder_index].unsqueeze(0))
            w_end = outputs_end[2]
            w_end = w_end.unsqueeze(axis=1)
            treat_predict_end = (outputs_end[0] * w_end) + (1 - w_end) * outputs_end[1]
            treat_predict_end = treat_predict_end.detach().cpu().numpy()
            temp = treat_predict_end[:, :window_size // 2].reshape(-1).tolist() + temp
            remainder_index = -1
            outputs_end = model(subseqOH[remainder_index].unsqueeze(0))
            w_end = outputs_end[2]
            w_end = w_end.unsqueeze(axis=1)
            treat_predict_end = (outputs_end[0] * w_end) + (1 - w_end) * outputs_end[1]
            treat_predict_end = treat_predict_end.detach().cpu().numpy()
            temp = treat_predict_end[:, :window_size // 2].reshape(-1).tolist() + temp

        transcriptPd = np.array(temp)
        print("Signal lenth of prediction and transcript length !", len(transcriptPd), len(seqTranscribed))
        if (len(transcriptPd) != len(seqTranscribed)):
            print("Error: signal lenth not equal transcript length !", len(transcriptPd), len(seqTranscribed))
            return False
    return transcriptPd


def predict1HeadFromSeq(seqTranscribed, window_size, model, mode, codeModel="OneHot", device=None):
    remainder = (len(seqTranscribed) % window_size)
    print(f"## remainder Length is: {remainder}")

    if mode == "Fragment":
        num_sublists = len(seqTranscribed) // window_size
        sublists = [seqTranscribed[i * window_size: (i + 1) * window_size] for i in range(num_sublists)]
        remainder = len(seqTranscribed) % window_size
        if remainder:
            sublists.append(seqTranscribed[-window_size:])

    elif mode == "midBase":
        sublists = [seqTranscribed[i: i + window_size] for i in np.arange(0, len(seqTranscribed) - window_size + 1, step=1)]

    if codeModel == "OneHot":
        vocab = list("AUGC")
        seqsInt = onehot_encode(sublists, vocab, 4)
        subseqOH = torch.from_numpy(np.asarray(seqsInt)).to(device)
    subseqOH = torch.from_numpy(np.asarray(seqsInt)).to(device)

    model.eval()
    with torch.no_grad():
        temp = []
        if (codeModel == "Embedding"):
            subseqOH = subseqOH.int()
        elif (codeModel == "OneHot"):
            subseqOH = subseqOH.float()

        if mode == "Fragment":
            if remainder:
                for subseqOHTurn in subseqOH[:-1].split(1000, 0):
                    treat_predict = model(subseqOHTurn)
                    treat_predict = treat_predict.detach().cpu().numpy().tolist()
                    temp.extend(treat_predict)
                temp = np.array(temp).reshape(-1).tolist()
                print("first temp len", len(temp))
                treat_predict_end = model(subseqOH[-1:])
                treat_predict_end = treat_predict_end.detach().cpu().numpy().tolist()
                print(remainder, len(temp), len(np.array(treat_predict_end).reshape(-1).tolist()[-remainder:]), len(np.array(treat_predict_end).reshape(-1).tolist()))
                temp = temp + np.array(treat_predict_end).reshape(-1).tolist()[-remainder:]
                print("remainder: ", remainder, len(treat_predict_end[-remainder:]), len(temp))
            else:
                for subseqOHTurn in subseqOH.split(1000, 0):
                    treat_predict = model(subseqOHTurn)
                    treat_predict = treat_predict.detach().cpu().numpy().tolist()
                    temp.extend(treat_predict)
                temp = np.array(temp).reshape(-1).tolist()

        elif mode == "midBase":
            for subseqOHTurn in subseqOH.split(1000, 0):
                treat_predict = model(subseqOHTurn)
                treat_predict = treat_predict[:, window_size // 2]
                treat_predict = treat_predict.detach().cpu().numpy().tolist()
                temp.extend(treat_predict)
            remainder_index = 0
            treat_predict_end = model(subseqOH[remainder_index].unsqueeze(0))
            treat_predict_end = treat_predict_end.detach().cpu().numpy()
            temp = treat_predict_end[:, :window_size // 2].reshape(-1).tolist() + temp
            remainder_index = -1
            treat_predict_end = model(subseqOH[remainder_index].unsqueeze(0))
            treat_predict_end = treat_predict_end.detach().cpu().numpy()
            temp = treat_predict_end[:, (window_size // 2 + 1):].reshape(-1).tolist() + temp

        transcriptPd = np.array(temp)
        print("Signal lenth of prediction and transcript length !", len(transcriptPd), len(seqTranscribed))
        if (len(transcriptPd) != len(seqTranscribed)):
            print("Error: signal lenth not equal transcript length !", len(transcriptPd), len(seqTranscribed))
            return False
    return transcriptPd
