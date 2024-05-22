;;;  -*- lexical-binding: t -*-
(require 'matrix)
(require 'ops)
(require 'seq)


(defun nn-layer (in out &optional activation)
  "Return a layer with weight matrix of shape IN x OUT and
bias with optional ACTIVATION function."
  (if activation 
      `(,(matrix-random in out) ,(matrix-random out 1) ,activation)
    `(,(matrix-random in out) ,(matrix-random out 1) identity)))


(defun nn--copy-bias-vector (bias times)
  "Copy the BIAS column vector in R^n TIMES number of times to
produce a matrix of shape n x TIMES, where each of the columns
 are the input BIAS column vector"
  (matrix-transpose (make-vector times (aref (matrix-transpose bias) 0))))

(defun nn--fill (shape value)
  "Return a matrix of SHAPE with all elements set to VALUE"
  (let ((nrows (aref shape 0))
        (ncols (aref shape 1)))
    (make-vector nrows (make-vector ncols value))))


(defun nn-forward (x layer)
  "Return the result of the forward pass on LAYER on input X"
  (let* ((input-shape (matrix-shape x))
         (batch-size (aref input-shape 1))
         (weight-matrix (nth 0 layer))
         (bias-matrix (nn--copy-bias-vector (nth 1 layer) batch-size))
         (activation (nth 2 layer)))
    
    (funcall activation (matrix-add (matrix-matmul (matrix-transpose weight-matrix) x) bias-matrix))))


(defun nn-forward-layers (x layers)
  "Return the result of the sequential forward pass through LAYERS
on input X"
  (seq-reduce #'nn-forward layers x))


(defun nn-crossentropy (label pred)
  "Calculate the crossentropy loss between LABEL and PRED"
  (let* ((sum (ops-reduce-sum (matrix-hadamard label (ops-log pred)) 1))
         (minus-one (nn--fill (matrix-shape sum) -1)))
    (matrix-hadamard minus-one sum)))

(defun nn--build-gradient (a-idx d-idx)
  (let ((batch-size (aref (matrix-shape d-idx) 0)))
    `(,(matrix-scalar-mul (/ 1.0 batch-size) (matrix-matmul a-idx d-idx))
      ,(matrix-transpose (ops-reduce-mean d-idx 1)))))

(defun nn-gradient (x y model)
  "Return the gradient of training example(s) X Y evaluated at MODEL
using crossentropy loss."
  (let* ((s (ops-softmax (nn-forward-layers x model))))
    (setq nn--DL (matrix-transpose (matrix-subtract s y))))

  (setq nn--D-matrices `(,nn--DL))
  (setq nn--A-matrices `(,x))

  ;; calculate the activation matrices
  (dotimes (idx (- (length model) 1))
    (setq nn--A-matrices (append nn--A-matrices `(,(nn-forward-layers x (seq-take model (+ 1 idx)))))))

  ;; calculate the d matrices
  (dotimes (count (- (length model) 1))
    (let* ((idx (- (- (length model) 1) count))
           (a-idx (nth idx nn--A-matrices))
           (h-idx (ops-heaviside (matrix-transpose a-idx)))
           (wt-idx-plus-one (matrix-transpose (nth 0 (nth idx model))))
           (d-idx-plus-one (car nn--D-matrices))
           (d-idx (matrix-hadamard (matrix-matmul d-idx-plus-one wt-idx-plus-one)
                                   h-idx)))
      (setq  nn--D-matrices (append `(,d-idx) nn--D-matrices))))


  (seq-mapn #'nn--build-gradient nn--A-matrices nn--D-matrices))


        

(provide 'nn)
