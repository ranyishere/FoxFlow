import torch
import numpy as np

# Seed for reproducibility
torch.manual_seed(42)

num_samples = 60000
images = torch.rand(num_samples, 1, 28, 28, dtype=torch.float64)
labels = torch.randint(0, 10, (num_samples,), dtype=torch.int64)

# Save as raw binary — compatible with LibTorch's torch::from_file
images.numpy().tofile("mnist_images.bin")
labels.numpy().tofile("mnist_labels.bin")

print(f"Saved images: {images.shape} dtype={images.dtype} -> mnist_images.bin")
print(f"Saved labels: {labels.shape} dtype={labels.dtype} -> mnist_labels.bin")

# Verify round-trip: load back and check shape matches what C++ will do
images_loaded = torch.from_numpy(np.fromfile("mnist_images.bin", dtype=np.float64)).reshape(60000, 1, 28, 28)
labels_loaded = torch.from_numpy(np.fromfile("mnist_labels.bin", dtype=np.int64)).reshape(60000)

batch_index = 0
print(f"\nIndex {batch_index} image_batch: {images_loaded[batch_index].shape}")
print(f"Index {batch_index} label:       {labels_loaded[batch_index].item()}")
print(f"\nBatch slice image_batch: {images_loaded[0:64].shape}")
print(f"Batch slice labels:      {labels_loaded[0:64].shape}")
