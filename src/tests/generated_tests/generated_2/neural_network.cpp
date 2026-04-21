#include <torch/torch.h>
#include <iostream>

// A simple 2-layer MLP: input -> hidden -> output
struct NetImpl : torch::nn::Module {
  torch::nn::Linear fc1{nullptr}, fc2{nullptr};

  NetImpl(int64_t in_dim, int64_t hidden_dim, int64_t out_dim) {
    fc1 = register_module("fc1", torch::nn::Linear(in_dim, hidden_dim));
    fc2 = register_module("fc2", torch::nn::Linear(hidden_dim, out_dim));
  }

  torch::Tensor forward(torch::Tensor x) {
    x = torch::relu(fc1->forward(x));
    x = fc2->forward(x);
    return x;
  }
};
TORCH_MODULE(Net); // creates Net as a shared_ptr wrapper around NetImpl

int main() {
  torch::manual_seed(0);

  // ----- 1) Create some fake regression data y = x0 - 2*x1 + noise
  const int64_t N = 2048;
  auto X = torch::randn({N, 2});
  auto y = (X.index({torch::indexing::Slice(), 0}) - 2 * X.index({torch::indexing::Slice(), 1}))
             .unsqueeze(1)
           + 0.1 * torch::randn({N, 1});

  // ----- 2) Model, loss, optimizer
  Net model(2, 64, 1);
  model->train();

  torch::optim::Adam optimizer(model->parameters(), torch::optim::AdamOptions(1e-3));
  auto loss_fn = torch::nn::MSELoss();

  // ----- 3) Training loop
  const int64_t epochs = 50;
  const int64_t batch_size = 64;

  for (int64_t epoch = 1; epoch <= epochs; ++epoch) {
    // shuffle indices each epoch
    auto idx = torch::randperm(N);
    double epoch_loss = 0.0;

    for (int64_t start = 0; start < N; start += batch_size) {
      auto end = std::min(start + batch_size, N);
      auto batch_idx = idx.slice(0, start, end);

      auto xb = X.index_select(0, batch_idx);
      auto yb = y.index_select(0, batch_idx);

      optimizer.zero_grad();

      auto pred = model->forward(xb);
      auto loss = loss_fn(pred, yb);

      loss.backward();
      optimizer.step();

      epoch_loss += loss.item<double>() * (end - start);
    }

    epoch_loss /= static_cast<double>(N);
    if (epoch % 10 == 0) {
      std::cout << "Epoch " << epoch << " | loss = " << epoch_loss << "\n";
    }
  }

  // ----- 4) Inference
  model->eval();
  auto test = torch::tensor({{1.0, 2.0}});
  auto out = model->forward(test);
  std::cout << "Pred( [1,2] ) = " << out << "\n";
}

