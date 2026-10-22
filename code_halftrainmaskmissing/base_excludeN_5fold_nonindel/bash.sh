python3 split_chr.py 0
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 0 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 1 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 2 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 3 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 4 
CUDA_VISIBLE_DEVICES=0 python3 predict.py 0 
CUDA_VISIBLE_DEVICES=0 python3 predict_train.py 0 

python3 split_chr.py 1 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 0 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 1 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 2 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 3 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 4 
CUDA_VISIBLE_DEVICES=0 python3 predict.py 1 
CUDA_VISIBLE_DEVICES=0 python3 predict_train.py 1 

python3 split_chr.py 2
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 0 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 1 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 2 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 3 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 4 
CUDA_VISIBLE_DEVICES=0 python3 predict.py 2 
CUDA_VISIBLE_DEVICES=0 python3 predict_train.py 2 

python3 split_chr.py 3
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 0 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 1 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 2 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 3 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 4 
CUDA_VISIBLE_DEVICES=0 python3 predict.py 3 
CUDA_VISIBLE_DEVICES=0 python3 predict_train.py 3 

python3 split_chr.py 4
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 0 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 1 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 2 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 3 
CUDA_VISIBLE_DEVICES=0 python3 train_batchnorm.py 4 
CUDA_VISIBLE_DEVICES=0 python3 predict.py 4 
CUDA_VISIBLE_DEVICES=0 python3 predict_train.py 4 


