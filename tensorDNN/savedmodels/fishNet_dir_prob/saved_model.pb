��/
��
B
AddV2
x"T
y"T
z"T"
Ttype:
2	��
B
AssignVariableOp
resource
value"dtype"
dtypetype�
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
�
Conv2D

input"T
filter"T
output"T"
Ttype:	
2"
strides	list(int)"
use_cudnn_on_gpubool(",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 "-
data_formatstringNHWC:
NHWCNCHW" 
	dilations	list(int)

.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
�
MaxPool

input"T
output"T"
Ttype0:
2	"
ksize	list(int)(0"
strides	list(int)(0",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 ":
data_formatstringNHWC:
NHWCNCHWNCHW_VECT_C
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(�
=
Mul
x"T
y"T
z"T"
Ttype:
2	�

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
E
Relu
features"T
activations"T"
Ttype:
2	
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
0
Sigmoid
x"T
y"T"
Ttype:

2
9
Softmax
logits"T
softmax"T"
Ttype:
2
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring �
@
StaticRegexFullMatch	
input

output
"
patternstring
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.4.32unknown8��*
�
conv2d_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameconv2d_3/kernel
{
#conv2d_3/kernel/Read/ReadVariableOpReadVariableOpconv2d_3/kernel*&
_output_shapes
:*
dtype0
r
conv2d_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv2d_3/bias
k
!conv2d_3/bias/Read/ReadVariableOpReadVariableOpconv2d_3/bias*
_output_shapes
:*
dtype0
�
conv2d_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameconv2d_4/kernel
{
#conv2d_4/kernel/Read/ReadVariableOpReadVariableOpconv2d_4/kernel*&
_output_shapes
:*
dtype0
r
conv2d_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv2d_4/bias
k
!conv2d_4/bias/Read/ReadVariableOpReadVariableOpconv2d_4/bias*
_output_shapes
:*
dtype0
�
conv2d_5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape: * 
shared_nameconv2d_5/kernel
{
#conv2d_5/kernel/Read/ReadVariableOpReadVariableOpconv2d_5/kernel*&
_output_shapes
: *
dtype0
r
conv2d_5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameconv2d_5/bias
k
!conv2d_5/bias/Read/ReadVariableOpReadVariableOpconv2d_5/bias*
_output_shapes
: *
dtype0
z
dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*
shared_namedense_3/kernel
s
"dense_3/kernel/Read/ReadVariableOpReadVariableOpdense_3/kernel* 
_output_shapes
:
��*
dtype0
q
dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*
shared_namedense_3/bias
j
 dense_3/bias/Read/ReadVariableOpReadVariableOpdense_3/bias*
_output_shapes	
:�*
dtype0
y
dense_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*
shared_namedense_4/kernel
r
"dense_4/kernel/Read/ReadVariableOpReadVariableOpdense_4/kernel*
_output_shapes
:	�*
dtype0
p
dense_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_4/bias
i
 dense_4/bias/Read/ReadVariableOpReadVariableOpdense_4/bias*
_output_shapes
:*
dtype0
x
dense_5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*
shared_namedense_5/kernel
q
"dense_5/kernel/Read/ReadVariableOpReadVariableOpdense_5/kernel*
_output_shapes

:*
dtype0
p
dense_5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_5/bias
i
 dense_5/bias/Read/ReadVariableOpReadVariableOpdense_5/bias*
_output_shapes
:*
dtype0
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0
h
VariableVarHandleOp*
_output_shapes
: *
dtype0	*
shape:*
shared_name
Variable
a
Variable/Read/ReadVariableOpReadVariableOpVariable*
_output_shapes
:*
dtype0	
l

Variable_1VarHandleOp*
_output_shapes
: *
dtype0	*
shape:*
shared_name
Variable_1
e
Variable_1/Read/ReadVariableOpReadVariableOp
Variable_1*
_output_shapes
:*
dtype0	
l

Variable_2VarHandleOp*
_output_shapes
: *
dtype0	*
shape:*
shared_name
Variable_2
e
Variable_2/Read/ReadVariableOpReadVariableOp
Variable_2*
_output_shapes
:*
dtype0	
l

Variable_3VarHandleOp*
_output_shapes
: *
dtype0	*
shape:*
shared_name
Variable_3
e
Variable_3/Read/ReadVariableOpReadVariableOp
Variable_3*
_output_shapes
:*
dtype0	
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
�
Adam/conv2d_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv2d_3/kernel/m
�
*Adam/conv2d_3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_3/kernel/m*&
_output_shapes
:*
dtype0
�
Adam/conv2d_3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv2d_3/bias/m
y
(Adam/conv2d_3/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_3/bias/m*
_output_shapes
:*
dtype0
�
Adam/conv2d_4/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv2d_4/kernel/m
�
*Adam/conv2d_4/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_4/kernel/m*&
_output_shapes
:*
dtype0
�
Adam/conv2d_4/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv2d_4/bias/m
y
(Adam/conv2d_4/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_4/bias/m*
_output_shapes
:*
dtype0
�
Adam/conv2d_5/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: *'
shared_nameAdam/conv2d_5/kernel/m
�
*Adam/conv2d_5/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_5/kernel/m*&
_output_shapes
: *
dtype0
�
Adam/conv2d_5/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/conv2d_5/bias/m
y
(Adam/conv2d_5/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_5/bias/m*
_output_shapes
: *
dtype0
�
Adam/dense_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*&
shared_nameAdam/dense_3/kernel/m
�
)Adam/dense_3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_3/kernel/m* 
_output_shapes
:
��*
dtype0

Adam/dense_3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*$
shared_nameAdam/dense_3/bias/m
x
'Adam/dense_3/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_3/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/dense_4/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*&
shared_nameAdam/dense_4/kernel/m
�
)Adam/dense_4/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_4/kernel/m*
_output_shapes
:	�*
dtype0
~
Adam/dense_4/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/dense_4/bias/m
w
'Adam/dense_4/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_4/bias/m*
_output_shapes
:*
dtype0
�
Adam/dense_5/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*&
shared_nameAdam/dense_5/kernel/m

)Adam/dense_5/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_5/kernel/m*
_output_shapes

:*
dtype0
~
Adam/dense_5/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/dense_5/bias/m
w
'Adam/dense_5/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_5/bias/m*
_output_shapes
:*
dtype0
�
Adam/conv2d_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv2d_3/kernel/v
�
*Adam/conv2d_3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_3/kernel/v*&
_output_shapes
:*
dtype0
�
Adam/conv2d_3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv2d_3/bias/v
y
(Adam/conv2d_3/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_3/bias/v*
_output_shapes
:*
dtype0
�
Adam/conv2d_4/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv2d_4/kernel/v
�
*Adam/conv2d_4/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_4/kernel/v*&
_output_shapes
:*
dtype0
�
Adam/conv2d_4/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv2d_4/bias/v
y
(Adam/conv2d_4/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_4/bias/v*
_output_shapes
:*
dtype0
�
Adam/conv2d_5/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: *'
shared_nameAdam/conv2d_5/kernel/v
�
*Adam/conv2d_5/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_5/kernel/v*&
_output_shapes
: *
dtype0
�
Adam/conv2d_5/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/conv2d_5/bias/v
y
(Adam/conv2d_5/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_5/bias/v*
_output_shapes
: *
dtype0
�
Adam/dense_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*&
shared_nameAdam/dense_3/kernel/v
�
)Adam/dense_3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_3/kernel/v* 
_output_shapes
:
��*
dtype0

Adam/dense_3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*$
shared_nameAdam/dense_3/bias/v
x
'Adam/dense_3/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_3/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/dense_4/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*&
shared_nameAdam/dense_4/kernel/v
�
)Adam/dense_4/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_4/kernel/v*
_output_shapes
:	�*
dtype0
~
Adam/dense_4/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/dense_4/bias/v
w
'Adam/dense_4/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_4/bias/v*
_output_shapes
:*
dtype0
�
Adam/dense_5/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*&
shared_nameAdam/dense_5/kernel/v

)Adam/dense_5/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_5/kernel/v*
_output_shapes

:*
dtype0
~
Adam/dense_5/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/dense_5/bias/v
w
'Adam/dense_5/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_5/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
�e
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�e
value�eB�e B�e
�
layer_with_weights-0
layer-0
layer-1
regularization_losses
trainable_variables
	variables
	keras_api

signatures
�
layer-0
	layer-1

layer_with_weights-0

layer-2
layer-3
layer-4
layer_with_weights-1
layer-5
layer-6
layer-7
layer_with_weights-2
layer-8
layer-9
layer-10
layer_with_weights-3
layer-11
layer-12
layer_with_weights-4
layer-13
layer-14
layer_with_weights-5
layer-15
	optimizer
regularization_losses
trainable_variables
	variables
	keras_api
R
regularization_losses
trainable_variables
	variables
 	keras_api
 
V
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
V
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
�
regularization_losses
-layer_metrics
.layer_regularization_losses
/metrics
trainable_variables

0layers
1non_trainable_variables
	variables
 
�
2layer-0
3layer-1
4layer-2
5layer-3
6regularization_losses
7trainable_variables
8	variables
9	keras_api

:	keras_api
h

!kernel
"bias
;regularization_losses
<trainable_variables
=	variables
>	keras_api
R
?regularization_losses
@trainable_variables
A	variables
B	keras_api
R
Cregularization_losses
Dtrainable_variables
E	variables
F	keras_api
h

#kernel
$bias
Gregularization_losses
Htrainable_variables
I	variables
J	keras_api
R
Kregularization_losses
Ltrainable_variables
M	variables
N	keras_api
R
Oregularization_losses
Ptrainable_variables
Q	variables
R	keras_api
h

%kernel
&bias
Sregularization_losses
Ttrainable_variables
U	variables
V	keras_api
R
Wregularization_losses
Xtrainable_variables
Y	variables
Z	keras_api
R
[regularization_losses
\trainable_variables
]	variables
^	keras_api
h

'kernel
(bias
_regularization_losses
`trainable_variables
a	variables
b	keras_api
R
cregularization_losses
dtrainable_variables
e	variables
f	keras_api
h

)kernel
*bias
gregularization_losses
htrainable_variables
i	variables
j	keras_api
R
kregularization_losses
ltrainable_variables
m	variables
n	keras_api
h

+kernel
,bias
oregularization_losses
ptrainable_variables
q	variables
r	keras_api
�
siter

tbeta_1

ubeta_2
	vdecay
wlearning_rate!m�"m�#m�$m�%m�&m�'m�(m�)m�*m�+m�,m�!v�"v�#v�$v�%v�&v�'v�(v�)v�*v�+v�,v�
 
V
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
V
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11
�
regularization_losses
xlayer_metrics
ylayer_regularization_losses
zmetrics
trainable_variables

{layers
|non_trainable_variables
	variables
 
 
 
�
regularization_losses
}layer_metrics
~layer_regularization_losses
metrics
trainable_variables
�layers
�non_trainable_variables
	variables
US
VARIABLE_VALUEconv2d_3/kernel0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEconv2d_3/bias0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv2d_4/kernel0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEconv2d_4/bias0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv2d_5/kernel0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEconv2d_5/bias0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUE
TR
VARIABLE_VALUEdense_3/kernel0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUE
RP
VARIABLE_VALUEdense_3/bias0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUE
TR
VARIABLE_VALUEdense_4/kernel0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUE
RP
VARIABLE_VALUEdense_4/bias0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEdense_5/kernel1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEdense_5/bias1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUE
 
 
 

0
1
 

	�_rng
�	keras_api

	�_rng
�	keras_api

	�_rng
�	keras_api

	�_rng
�	keras_api
 
 
 
�
6regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
7trainable_variables
�layers
�non_trainable_variables
8	variables
 
 

!0
"1

!0
"1
�
;regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
<trainable_variables
�layers
�non_trainable_variables
=	variables
 
 
 
�
?regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
@trainable_variables
�layers
�non_trainable_variables
A	variables
 
 
 
�
Cregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Dtrainable_variables
�layers
�non_trainable_variables
E	variables
 

#0
$1

#0
$1
�
Gregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Htrainable_variables
�layers
�non_trainable_variables
I	variables
 
 
 
�
Kregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Ltrainable_variables
�layers
�non_trainable_variables
M	variables
 
 
 
�
Oregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Ptrainable_variables
�layers
�non_trainable_variables
Q	variables
 

%0
&1

%0
&1
�
Sregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Ttrainable_variables
�layers
�non_trainable_variables
U	variables
 
 
 
�
Wregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Xtrainable_variables
�layers
�non_trainable_variables
Y	variables
 
 
 
�
[regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
\trainable_variables
�layers
�non_trainable_variables
]	variables
 

'0
(1

'0
(1
�
_regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
`trainable_variables
�layers
�non_trainable_variables
a	variables
 
 
 
�
cregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
dtrainable_variables
�layers
�non_trainable_variables
e	variables
 

)0
*1

)0
*1
�
gregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
htrainable_variables
�layers
�non_trainable_variables
i	variables
 
 
 
�
kregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
ltrainable_variables
�layers
�non_trainable_variables
m	variables
 

+0
,1

+0
,1
�
oregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
ptrainable_variables
�layers
�non_trainable_variables
q	variables
][
VARIABLE_VALUE	Adam/iter>layer_with_weights-0/optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
a_
VARIABLE_VALUEAdam/beta_1@layer_with_weights-0/optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
a_
VARIABLE_VALUEAdam/beta_2@layer_with_weights-0/optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
_]
VARIABLE_VALUE
Adam/decay?layer_with_weights-0/optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUEAdam/learning_rateGlayer_with_weights-0/optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
 
 

�0
�1
v
0
	1

2
3
4
5
6
7
8
9
10
11
12
13
14
15
 
 
 
 
 
 

�
_state_var
 

�
_state_var
 

�
_state_var
 

�
_state_var
 
 
 
 

20
31
42
53
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
8

�total

�count
�	variables
�	keras_api
I

�total

�count
�
_fn_kwargs
�	variables
�	keras_api
mk
VARIABLE_VALUEVariableOlayer_with_weights-0/layer-0/layer-0/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUE
Variable_1Olayer_with_weights-0/layer-0/layer-1/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUE
Variable_2Olayer_with_weights-0/layer-0/layer-2/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUE
om
VARIABLE_VALUE
Variable_3Olayer_with_weights-0/layer-0/layer-3/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUE
db
VARIABLE_VALUEtotalIlayer_with_weights-0/keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
db
VARIABLE_VALUEcountIlayer_with_weights-0/keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

�0
�1

�	variables
fd
VARIABLE_VALUEtotal_1Ilayer_with_weights-0/keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
fd
VARIABLE_VALUEcount_1Ilayer_with_weights-0/keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE
 

�0
�1

�	variables
��
VARIABLE_VALUEAdam/conv2d_3/kernel/matrainable_variables/0/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_3/bias/matrainable_variables/1/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_4/kernel/matrainable_variables/2/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_4/bias/matrainable_variables/3/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_5/kernel/matrainable_variables/4/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_5/bias/matrainable_variables/5/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_3/kernel/matrainable_variables/6/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_3/bias/matrainable_variables/7/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_4/kernel/matrainable_variables/8/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_4/bias/matrainable_variables/9/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_5/kernel/mbtrainable_variables/10/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_5/bias/mbtrainable_variables/11/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_3/kernel/vatrainable_variables/0/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_3/bias/vatrainable_variables/1/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_4/kernel/vatrainable_variables/2/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_4/bias/vatrainable_variables/3/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_5/kernel/vatrainable_variables/4/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/conv2d_5/bias/vatrainable_variables/5/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_3/kernel/vatrainable_variables/6/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_3/bias/vatrainable_variables/7/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_4/kernel/vatrainable_variables/8/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_4/bias/vatrainable_variables/9/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_5/kernel/vbtrainable_variables/10/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
��
VARIABLE_VALUEAdam/dense_5/bias/vbtrainable_variables/11/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
�
"serving_default_sequential_4_inputPlaceholder*/
_output_shapes
:���������&*
dtype0*$
shape:���������&
�
StatefulPartitionedCallStatefulPartitionedCall"serving_default_sequential_4_inputconv2d_3/kernelconv2d_3/biasconv2d_4/kernelconv2d_4/biasconv2d_5/kernelconv2d_5/biasdense_3/kerneldense_3/biasdense_4/kerneldense_4/biasdense_5/kerneldense_5/bias*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *,
f'R%
#__inference_signature_wrapper_32315
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename#conv2d_3/kernel/Read/ReadVariableOp!conv2d_3/bias/Read/ReadVariableOp#conv2d_4/kernel/Read/ReadVariableOp!conv2d_4/bias/Read/ReadVariableOp#conv2d_5/kernel/Read/ReadVariableOp!conv2d_5/bias/Read/ReadVariableOp"dense_3/kernel/Read/ReadVariableOp dense_3/bias/Read/ReadVariableOp"dense_4/kernel/Read/ReadVariableOp dense_4/bias/Read/ReadVariableOp"dense_5/kernel/Read/ReadVariableOp dense_5/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOpVariable/Read/ReadVariableOpVariable_1/Read/ReadVariableOpVariable_2/Read/ReadVariableOpVariable_3/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOp*Adam/conv2d_3/kernel/m/Read/ReadVariableOp(Adam/conv2d_3/bias/m/Read/ReadVariableOp*Adam/conv2d_4/kernel/m/Read/ReadVariableOp(Adam/conv2d_4/bias/m/Read/ReadVariableOp*Adam/conv2d_5/kernel/m/Read/ReadVariableOp(Adam/conv2d_5/bias/m/Read/ReadVariableOp)Adam/dense_3/kernel/m/Read/ReadVariableOp'Adam/dense_3/bias/m/Read/ReadVariableOp)Adam/dense_4/kernel/m/Read/ReadVariableOp'Adam/dense_4/bias/m/Read/ReadVariableOp)Adam/dense_5/kernel/m/Read/ReadVariableOp'Adam/dense_5/bias/m/Read/ReadVariableOp*Adam/conv2d_3/kernel/v/Read/ReadVariableOp(Adam/conv2d_3/bias/v/Read/ReadVariableOp*Adam/conv2d_4/kernel/v/Read/ReadVariableOp(Adam/conv2d_4/bias/v/Read/ReadVariableOp*Adam/conv2d_5/kernel/v/Read/ReadVariableOp(Adam/conv2d_5/bias/v/Read/ReadVariableOp)Adam/dense_3/kernel/v/Read/ReadVariableOp'Adam/dense_3/bias/v/Read/ReadVariableOp)Adam/dense_4/kernel/v/Read/ReadVariableOp'Adam/dense_4/bias/v/Read/ReadVariableOp)Adam/dense_5/kernel/v/Read/ReadVariableOp'Adam/dense_5/bias/v/Read/ReadVariableOpConst*>
Tin7
523					*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *'
f"R 
__inference__traced_save_34286
�	
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameconv2d_3/kernelconv2d_3/biasconv2d_4/kernelconv2d_4/biasconv2d_5/kernelconv2d_5/biasdense_3/kerneldense_3/biasdense_4/kerneldense_4/biasdense_5/kerneldense_5/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_rateVariable
Variable_1
Variable_2
Variable_3totalcounttotal_1count_1Adam/conv2d_3/kernel/mAdam/conv2d_3/bias/mAdam/conv2d_4/kernel/mAdam/conv2d_4/bias/mAdam/conv2d_5/kernel/mAdam/conv2d_5/bias/mAdam/dense_3/kernel/mAdam/dense_3/bias/mAdam/dense_4/kernel/mAdam/dense_4/bias/mAdam/dense_5/kernel/mAdam/dense_5/bias/mAdam/conv2d_3/kernel/vAdam/conv2d_3/bias/vAdam/conv2d_4/kernel/vAdam/conv2d_4/bias/vAdam/conv2d_5/kernel/vAdam/conv2d_5/bias/vAdam/dense_3/kernel/vAdam/dense_3/bias/vAdam/dense_4/kernel/vAdam/dense_4/bias/vAdam/dense_5/kernel/vAdam/dense_5/bias/v*=
Tin6
422*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� **
f%R#
!__inference__traced_restore_34443�(
�
b
D__inference_dropout_6_layer_call_and_return_conditional_losses_31586

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:����������2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:����������2

Identity_1"!

identity_1Identity_1:output:0*'
_input_shapes
:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�

�
C__inference_conv2d_4_layer_call_and_return_conditional_losses_31447

inputs"
conv2d_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������*
paddingSAME*
strides
2
Conv2D�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*
T0*/
_output_shapes
:���������2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
�
E
)__inference_dropout_5_layer_call_fn_33905

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_5_layer_call_and_return_conditional_losses_314802
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
�
b
)__inference_dropout_4_layer_call_fn_33853

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_4_layer_call_and_return_conditional_losses_314172
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
E
)__inference_softmax_1_layer_call_fn_33463

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_softmax_1_layer_call_and_return_conditional_losses_320542
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*&
_input_shapes
:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
b
D__inference_dropout_4_layer_call_and_return_conditional_losses_33848

inputs

identity_1b
IdentityIdentityinputs*
T0*/
_output_shapes
:���������&2

Identityq

Identity_1IdentityIdentity:output:0*
T0*/
_output_shapes
:���������&2

Identity_1"!

identity_1Identity_1:output:0*.
_input_shapes
:���������&:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
E
)__inference_dropout_7_layer_call_fn_34054

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_7_layer_call_and_return_conditional_losses_316492
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*&
_input_shapes
:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
B__inference_dense_4_layer_call_and_return_conditional_losses_31616

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
E
)__inference_flatten_1_layer_call_fn_33936

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_flatten_1_layer_call_and_return_conditional_losses_315282
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*.
_input_shapes
:��������� :W S
/
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
K
/__inference_max_pooling2d_4_layer_call_fn_31335

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4������������������������������������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_4_layer_call_and_return_conditional_losses_313292
PartitionedCall�
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4������������������������������������2

Identity"
identityIdentity:output:0*I
_input_shapes8
6:4������������������������������������:r n
J
_output_shapes8
6:4������������������������������������
 
_user_specified_nameinputs
�
�
B__inference_dense_3_layer_call_and_return_conditional_losses_31553

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�0dense_3/kernel/Regularizer/Square/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:����������2
Relu�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�	
�
,__inference_sequential_5_layer_call_fn_32875

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_5_layer_call_and_return_conditional_losses_322452
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
|
'__inference_dense_3_layer_call_fn_33968

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_315532
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*/
_input_shapes
:����������::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
c
D__inference_dropout_5_layer_call_and_return_conditional_losses_33890

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Const{
dropout/MulMulinputsdropout/Const:output:0*
T0*/
_output_shapes
:���������2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*/
_output_shapes
:���������*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*/
_output_shapes
:���������2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*/
_output_shapes
:���������2
dropout/Cast�
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*/
_output_shapes
:���������2
dropout/Mul_1m
IdentityIdentitydropout/Mul_1:z:0*
T0*/
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
,__inference_sequential_5_layer_call_fn_32201
sequential_4_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallsequential_4_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_5_layer_call_and_return_conditional_losses_321682
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_4_input
�"
�
G__inference_sequential_5_layer_call_and_return_conditional_losses_32245

inputs
sequential_4_32206
sequential_4_32208
sequential_4_32210
sequential_4_32212
sequential_4_32214
sequential_4_32216
sequential_4_32218
sequential_4_32220
sequential_4_32222
sequential_4_32224
sequential_4_32226
sequential_4_32228
identity��0dense_3/kernel/Regularizer/Square/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�$sequential_4/StatefulPartitionedCall�
$sequential_4/StatefulPartitionedCallStatefulPartitionedCallinputssequential_4_32206sequential_4_32208sequential_4_32210sequential_4_32212sequential_4_32214sequential_4_32216sequential_4_32218sequential_4_32220sequential_4_32222sequential_4_32224sequential_4_32226sequential_4_32228*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_4_layer_call_and_return_conditional_losses_319232&
$sequential_4/StatefulPartitionedCall�
softmax_1/PartitionedCallPartitionedCall-sequential_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_softmax_1_layer_call_and_return_conditional_losses_320542
softmax_1/PartitionedCall�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_32218* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_32222*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity"softmax_1/PartitionedCall:output:01^dense_3/kernel/Regularizer/Square/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp%^sequential_4/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
c
D__inference_dropout_7_layer_call_and_return_conditional_losses_31644

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:���������2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:���������*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���>2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:���������2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:���������2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:���������2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*&
_input_shapes
:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
c
D__inference_dropout_4_layer_call_and_return_conditional_losses_33843

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Const{
dropout/MulMulinputsdropout/Const:output:0*
T0*/
_output_shapes
:���������&2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*/
_output_shapes
:���������&*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*/
_output_shapes
:���������&2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*/
_output_shapes
:���������&2
dropout/Cast�
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*/
_output_shapes
:���������&2
dropout/Mul_1m
IdentityIdentitydropout/Mul_1:z:0*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
E
)__inference_dropout_4_layer_call_fn_33858

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_4_layer_call_and_return_conditional_losses_314222
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
b
)__inference_dropout_7_layer_call_fn_34049

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_7_layer_call_and_return_conditional_losses_316442
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*&
_input_shapes
:���������22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
C__inference_conv2d_4_layer_call_and_return_conditional_losses_33869

inputs"
conv2d_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������*
paddingSAME*
strides
2
Conv2D�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:���������2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*
T0*/
_output_shapes
:���������2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
B__inference_dense_5_layer_call_and_return_conditional_losses_31673

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAdda
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:���������2	
Sigmoid�
IdentityIdentitySigmoid:y:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
H
,__inference_sequential_3_layer_call_fn_33811

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_313082
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
`
D__inference_softmax_1_layer_call_and_return_conditional_losses_32054

inputs
identityW
SoftmaxSoftmaxinputs*
T0*'
_output_shapes
:���������2	
Softmaxe
IdentityIdentitySoftmax:softmax:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*&
_input_shapes
:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
C__inference_conv2d_5_layer_call_and_return_conditional_losses_33916

inputs"
conv2d_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: *
dtype02
Conv2D/ReadVariableOp�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 *
paddingSAME*
strides
2
Conv2D�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:���������	 2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*
T0*/
_output_shapes
:���������	 2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������	::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:���������	
 
_user_specified_nameinputs
�
b
D__inference_dropout_7_layer_call_and_return_conditional_losses_34044

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:���������2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:���������2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
b
D__inference_dropout_5_layer_call_and_return_conditional_losses_31480

inputs

identity_1b
IdentityIdentityinputs*
T0*/
_output_shapes
:���������2

Identityq

Identity_1IdentityIdentity:output:0*
T0*/
_output_shapes
:���������2

Identity_1"!

identity_1Identity_1:output:0*.
_input_shapes
:���������:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
��
�
G__inference_sequential_5_layer_call_and_return_conditional_losses_32739

inputsY
Usequential_4_sequential_3_random_rotation_1_stateful_uniform_statefuluniform_resource\
Xsequential_4_sequential_3_random_translation_1_stateful_uniform_statefuluniform_resourceU
Qsequential_4_sequential_3_random_zoom_1_stateful_uniform_statefuluniform_resource8
4sequential_4_conv2d_3_conv2d_readvariableop_resource9
5sequential_4_conv2d_3_biasadd_readvariableop_resource8
4sequential_4_conv2d_4_conv2d_readvariableop_resource9
5sequential_4_conv2d_4_biasadd_readvariableop_resource8
4sequential_4_conv2d_5_conv2d_readvariableop_resource9
5sequential_4_conv2d_5_biasadd_readvariableop_resource7
3sequential_4_dense_3_matmul_readvariableop_resource8
4sequential_4_dense_3_biasadd_readvariableop_resource7
3sequential_4_dense_4_matmul_readvariableop_resource8
4sequential_4_dense_4_biasadd_readvariableop_resource7
3sequential_4_dense_5_matmul_readvariableop_resource8
4sequential_4_dense_5_biasadd_readvariableop_resource
identity��0dense_3/kernel/Regularizer/Square/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�,sequential_4/conv2d_3/BiasAdd/ReadVariableOp�+sequential_4/conv2d_3/Conv2D/ReadVariableOp�,sequential_4/conv2d_4/BiasAdd/ReadVariableOp�+sequential_4/conv2d_4/Conv2D/ReadVariableOp�,sequential_4/conv2d_5/BiasAdd/ReadVariableOp�+sequential_4/conv2d_5/Conv2D/ReadVariableOp�+sequential_4/dense_3/BiasAdd/ReadVariableOp�*sequential_4/dense_3/MatMul/ReadVariableOp�+sequential_4/dense_4/BiasAdd/ReadVariableOp�*sequential_4/dense_4/MatMul/ReadVariableOp�+sequential_4/dense_5/BiasAdd/ReadVariableOp�*sequential_4/dense_5/MatMul/ReadVariableOp�Lsequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniform�Osequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniform�Qsequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniform�Hsequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniform�
Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/control_dependencyIdentityinputs*
T0*
_class
loc:@inputs*/
_output_shapes
:���������&2S
Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/control_dependency�
Dsequential_4/sequential_3/random_flip_1/random_flip_left_right/ShapeShapeZsequential_4/sequential_3/random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:2F
Dsequential_4/sequential_3/random_flip_1/random_flip_left_right/Shape�
Rsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2T
Rsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stack�
Tsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2V
Tsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_1�
Tsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2V
Tsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_2�
Lsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_sliceStridedSliceMsequential_4/sequential_3/random_flip_1/random_flip_left_right/Shape:output:0[sequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stack:output:0]sequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_1:output:0]sequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2N
Lsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice�
Ssequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/shapePackUsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2U
Ssequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/shape�
Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2S
Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/min�
Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2S
Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/max�
[sequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/RandomUniformRandomUniform\sequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/shape:output:0*
T0*#
_output_shapes
:���������*
dtype02]
[sequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/RandomUniform�
Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/MulMuldsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/RandomUniform:output:0Zsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/max:output:0*
T0*#
_output_shapes
:���������2S
Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/Mul�
Nsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2P
Nsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/1�
Nsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2P
Nsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/2�
Nsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2P
Nsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/3�
Lsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shapePackUsequential_4/sequential_3/random_flip_1/random_flip_left_right/strided_slice:output:0Wsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/1:output:0Wsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/2:output:0Wsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2N
Lsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape�
Fsequential_4/sequential_3/random_flip_1/random_flip_left_right/ReshapeReshapeUsequential_4/sequential_3/random_flip_1/random_flip_left_right/random_uniform/Mul:z:0Usequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:���������2H
Fsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape�
Dsequential_4/sequential_3/random_flip_1/random_flip_left_right/RoundRoundOsequential_4/sequential_3/random_flip_1/random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:���������2F
Dsequential_4/sequential_3/random_flip_1/random_flip_left_right/Round�
Msequential_4/sequential_3/random_flip_1/random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:2O
Msequential_4/sequential_3/random_flip_1/random_flip_left_right/ReverseV2/axis�
Hsequential_4/sequential_3/random_flip_1/random_flip_left_right/ReverseV2	ReverseV2Zsequential_4/sequential_3/random_flip_1/random_flip_left_right/control_dependency:output:0Vsequential_4/sequential_3/random_flip_1/random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:���������&2J
Hsequential_4/sequential_3/random_flip_1/random_flip_left_right/ReverseV2�
Bsequential_4/sequential_3/random_flip_1/random_flip_left_right/mulMulHsequential_4/sequential_3/random_flip_1/random_flip_left_right/Round:y:0Qsequential_4/sequential_3/random_flip_1/random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:���������&2D
Bsequential_4/sequential_3/random_flip_1/random_flip_left_right/mul�
Dsequential_4/sequential_3/random_flip_1/random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2F
Dsequential_4/sequential_3/random_flip_1/random_flip_left_right/sub/x�
Bsequential_4/sequential_3/random_flip_1/random_flip_left_right/subSubMsequential_4/sequential_3/random_flip_1/random_flip_left_right/sub/x:output:0Hsequential_4/sequential_3/random_flip_1/random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:���������2D
Bsequential_4/sequential_3/random_flip_1/random_flip_left_right/sub�
Dsequential_4/sequential_3/random_flip_1/random_flip_left_right/mul_1MulFsequential_4/sequential_3/random_flip_1/random_flip_left_right/sub:z:0Zsequential_4/sequential_3/random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:���������&2F
Dsequential_4/sequential_3/random_flip_1/random_flip_left_right/mul_1�
Bsequential_4/sequential_3/random_flip_1/random_flip_left_right/addAddV2Fsequential_4/sequential_3/random_flip_1/random_flip_left_right/mul:z:0Hsequential_4/sequential_3/random_flip_1/random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:���������&2D
Bsequential_4/sequential_3/random_flip_1/random_flip_left_right/add�
1sequential_4/sequential_3/random_rotation_1/ShapeShapeFsequential_4/sequential_3/random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:23
1sequential_4/sequential_3/random_rotation_1/Shape�
?sequential_4/sequential_3/random_rotation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2A
?sequential_4/sequential_3/random_rotation_1/strided_slice/stack�
Asequential_4/sequential_3/random_rotation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2C
Asequential_4/sequential_3/random_rotation_1/strided_slice/stack_1�
Asequential_4/sequential_3/random_rotation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2C
Asequential_4/sequential_3/random_rotation_1/strided_slice/stack_2�
9sequential_4/sequential_3/random_rotation_1/strided_sliceStridedSlice:sequential_4/sequential_3/random_rotation_1/Shape:output:0Hsequential_4/sequential_3/random_rotation_1/strided_slice/stack:output:0Jsequential_4/sequential_3/random_rotation_1/strided_slice/stack_1:output:0Jsequential_4/sequential_3/random_rotation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2;
9sequential_4/sequential_3/random_rotation_1/strided_slice�
Asequential_4/sequential_3/random_rotation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2C
Asequential_4/sequential_3/random_rotation_1/strided_slice_1/stack�
Csequential_4/sequential_3/random_rotation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2E
Csequential_4/sequential_3/random_rotation_1/strided_slice_1/stack_1�
Csequential_4/sequential_3/random_rotation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2E
Csequential_4/sequential_3/random_rotation_1/strided_slice_1/stack_2�
;sequential_4/sequential_3/random_rotation_1/strided_slice_1StridedSlice:sequential_4/sequential_3/random_rotation_1/Shape:output:0Jsequential_4/sequential_3/random_rotation_1/strided_slice_1/stack:output:0Lsequential_4/sequential_3/random_rotation_1/strided_slice_1/stack_1:output:0Lsequential_4/sequential_3/random_rotation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2=
;sequential_4/sequential_3/random_rotation_1/strided_slice_1�
0sequential_4/sequential_3/random_rotation_1/CastCastDsequential_4/sequential_3/random_rotation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 22
0sequential_4/sequential_3/random_rotation_1/Cast�
Asequential_4/sequential_3/random_rotation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2C
Asequential_4/sequential_3/random_rotation_1/strided_slice_2/stack�
Csequential_4/sequential_3/random_rotation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2E
Csequential_4/sequential_3/random_rotation_1/strided_slice_2/stack_1�
Csequential_4/sequential_3/random_rotation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2E
Csequential_4/sequential_3/random_rotation_1/strided_slice_2/stack_2�
;sequential_4/sequential_3/random_rotation_1/strided_slice_2StridedSlice:sequential_4/sequential_3/random_rotation_1/Shape:output:0Jsequential_4/sequential_3/random_rotation_1/strided_slice_2/stack:output:0Lsequential_4/sequential_3/random_rotation_1/strided_slice_2/stack_1:output:0Lsequential_4/sequential_3/random_rotation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2=
;sequential_4/sequential_3/random_rotation_1/strided_slice_2�
2sequential_4/sequential_3/random_rotation_1/Cast_1CastDsequential_4/sequential_3/random_rotation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 24
2sequential_4/sequential_3/random_rotation_1/Cast_1�
Bsequential_4/sequential_3/random_rotation_1/stateful_uniform/shapePackBsequential_4/sequential_3/random_rotation_1/strided_slice:output:0*
N*
T0*
_output_shapes
:2D
Bsequential_4/sequential_3/random_rotation_1/stateful_uniform/shape�
@sequential_4/sequential_3/random_rotation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *   �2B
@sequential_4/sequential_3/random_rotation_1/stateful_uniform/min�
@sequential_4/sequential_3/random_rotation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2B
@sequential_4/sequential_3/random_rotation_1/stateful_uniform/max�
Vsequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2X
Vsequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniform/algorithm�
Lsequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniformStatefulUniformUsequential_4_sequential_3_random_rotation_1_stateful_uniform_statefuluniform_resource_sequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniform/algorithm:output:0Ksequential_4/sequential_3/random_rotation_1/stateful_uniform/shape:output:0*#
_output_shapes
:���������*
shape_dtype02N
Lsequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniform�
@sequential_4/sequential_3/random_rotation_1/stateful_uniform/subSubIsequential_4/sequential_3/random_rotation_1/stateful_uniform/max:output:0Isequential_4/sequential_3/random_rotation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2B
@sequential_4/sequential_3/random_rotation_1/stateful_uniform/sub�
@sequential_4/sequential_3/random_rotation_1/stateful_uniform/mulMulUsequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniform:output:0Dsequential_4/sequential_3/random_rotation_1/stateful_uniform/sub:z:0*
T0*#
_output_shapes
:���������2B
@sequential_4/sequential_3/random_rotation_1/stateful_uniform/mul�
<sequential_4/sequential_3/random_rotation_1/stateful_uniformAddDsequential_4/sequential_3/random_rotation_1/stateful_uniform/mul:z:0Isequential_4/sequential_3/random_rotation_1/stateful_uniform/min:output:0*
T0*#
_output_shapes
:���������2>
<sequential_4/sequential_3/random_rotation_1/stateful_uniform�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub/y�
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/subSub6sequential_4/sequential_3/random_rotation_1/Cast_1:y:0Jsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 2A
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/sub�
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/CosCos@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2A
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos�
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2E
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_1/y�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_1Sub6sequential_4/sequential_3/random_rotation_1/Cast_1:y:0Lsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_1�
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/mulMulCsequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos:y:0Esequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:���������2A
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/mul�
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/SinSin@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2A
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin�
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2E
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_2/y�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_2Sub4sequential_4/sequential_3/random_rotation_1/Cast:y:0Lsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_2�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_1MulCsequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin:y:0Esequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_1�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_3SubCsequential_4/sequential_3/random_rotation_1/rotation_matrix/mul:z:0Esequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_3�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_4SubCsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub:z:0Esequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_4�
Esequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2G
Esequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv/y�
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/truedivRealDivEsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_4:z:0Nsequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:���������2E
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv�
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2E
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_5/y�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_5Sub4sequential_4/sequential_3/random_rotation_1/Cast:y:0Lsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_5�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_1Sin@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_1�
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2E
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_6/y�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_6Sub6sequential_4/sequential_3/random_rotation_1/Cast_1:y:0Lsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_6�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_2MulEsequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_1:y:0Esequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_2�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_1Cos@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_1�
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2E
Csequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_7/y�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_7Sub4sequential_4/sequential_3/random_rotation_1/Cast:y:0Lsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_7�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_3MulEsequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_1:y:0Esequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_3�
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/addAddV2Esequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_2:z:0Esequential_4/sequential_3/random_rotation_1/rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:���������2A
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/add�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_8SubEsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_5:z:0Csequential_4/sequential_3/random_rotation_1/rotation_matrix/add:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_8�
Gsequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2I
Gsequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv_1/y�
Esequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv_1RealDivEsequential_4/sequential_3/random_rotation_1/rotation_matrix/sub_8:z:0Psequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:���������2G
Esequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv_1�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/ShapeShape@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*
_output_shapes
:2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Shape�
Osequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2Q
Osequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stack�
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2S
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_1�
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2S
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_2�
Isequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_sliceStridedSliceJsequential_4/sequential_3/random_rotation_1/rotation_matrix/Shape:output:0Xsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stack:output:0Zsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_1:output:0Zsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2K
Isequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_2Cos@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_2�
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        2S
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_1�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_2�
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1StridedSliceEsequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_2:y:0Zsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_1:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2M
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_2Sin@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_2�
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        2S
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_1�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_2�
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2StridedSliceEsequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_2:y:0Zsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_1:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2M
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2�
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/NegNegTsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2A
?sequential_4/sequential_3/random_rotation_1/rotation_matrix/Neg�
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        2S
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_1�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_2�
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3StridedSliceGsequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv:z:0Zsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_1:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2M
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_3Sin@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_3�
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        2S
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_1�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_2�
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4StridedSliceEsequential_4/sequential_3/random_rotation_1/rotation_matrix/Sin_3:y:0Zsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_1:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2M
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_3Cos@sequential_4/sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_3�
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        2S
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_1�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_2�
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5StridedSliceEsequential_4/sequential_3/random_rotation_1/rotation_matrix/Cos_3:y:0Zsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_1:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2M
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5�
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        2S
Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_1�
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2U
Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_2�
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6StridedSliceIsequential_4/sequential_3/random_rotation_1/rotation_matrix/truediv_1:z:0Zsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_1:output:0\sequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2M
Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6�
Gsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2I
Gsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/mul/y�
Esequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/mulMulRsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice:output:0Psequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2G
Esequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/mul�
Hsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2J
Hsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/Less/y�
Fsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/LessLessIsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/mul:z:0Qsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2H
Fsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/Less�
Jsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2L
Jsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/packed/1�
Hsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/packedPackRsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice:output:0Ssequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2J
Hsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/packed�
Gsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2I
Gsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/Const�
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/zerosFillQsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/packed:output:0Psequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2C
Asequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros�
Gsequential_4/sequential_3/random_rotation_1/rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2I
Gsequential_4/sequential_3/random_rotation_1/rotation_matrix/concat/axis�
Bsequential_4/sequential_3/random_rotation_1/rotation_matrix/concatConcatV2Tsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_1:output:0Csequential_4/sequential_3/random_rotation_1/rotation_matrix/Neg:y:0Tsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_3:output:0Tsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_4:output:0Tsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_5:output:0Tsequential_4/sequential_3/random_rotation_1/rotation_matrix/strided_slice_6:output:0Jsequential_4/sequential_3/random_rotation_1/rotation_matrix/zeros:output:0Psequential_4/sequential_3/random_rotation_1/rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2D
Bsequential_4/sequential_3/random_rotation_1/rotation_matrix/concat�
;sequential_4/sequential_3/random_rotation_1/transform/ShapeShapeFsequential_4/sequential_3/random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:2=
;sequential_4/sequential_3/random_rotation_1/transform/Shape�
Isequential_4/sequential_3/random_rotation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2K
Isequential_4/sequential_3/random_rotation_1/transform/strided_slice/stack�
Ksequential_4/sequential_3/random_rotation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2M
Ksequential_4/sequential_3/random_rotation_1/transform/strided_slice/stack_1�
Ksequential_4/sequential_3/random_rotation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2M
Ksequential_4/sequential_3/random_rotation_1/transform/strided_slice/stack_2�
Csequential_4/sequential_3/random_rotation_1/transform/strided_sliceStridedSliceDsequential_4/sequential_3/random_rotation_1/transform/Shape:output:0Rsequential_4/sequential_3/random_rotation_1/transform/strided_slice/stack:output:0Tsequential_4/sequential_3/random_rotation_1/transform/strided_slice/stack_1:output:0Tsequential_4/sequential_3/random_rotation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2E
Csequential_4/sequential_3/random_rotation_1/transform/strided_slice�
@sequential_4/sequential_3/random_rotation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2B
@sequential_4/sequential_3/random_rotation_1/transform/fill_value�
Psequential_4/sequential_3/random_rotation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Fsequential_4/sequential_3/random_flip_1/random_flip_left_right/add:z:0Ksequential_4/sequential_3/random_rotation_1/rotation_matrix/concat:output:0Lsequential_4/sequential_3/random_rotation_1/transform/strided_slice:output:0Isequential_4/sequential_3/random_rotation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2R
Psequential_4/sequential_3/random_rotation_1/transform/ImageProjectiveTransformV3�
4sequential_4/sequential_3/random_translation_1/ShapeShapeesequential_4/sequential_3/random_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:26
4sequential_4/sequential_3/random_translation_1/Shape�
Bsequential_4/sequential_3/random_translation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2D
Bsequential_4/sequential_3/random_translation_1/strided_slice/stack�
Dsequential_4/sequential_3/random_translation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2F
Dsequential_4/sequential_3/random_translation_1/strided_slice/stack_1�
Dsequential_4/sequential_3/random_translation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2F
Dsequential_4/sequential_3/random_translation_1/strided_slice/stack_2�
<sequential_4/sequential_3/random_translation_1/strided_sliceStridedSlice=sequential_4/sequential_3/random_translation_1/Shape:output:0Ksequential_4/sequential_3/random_translation_1/strided_slice/stack:output:0Msequential_4/sequential_3/random_translation_1/strided_slice/stack_1:output:0Msequential_4/sequential_3/random_translation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2>
<sequential_4/sequential_3/random_translation_1/strided_slice�
Dsequential_4/sequential_3/random_translation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2F
Dsequential_4/sequential_3/random_translation_1/strided_slice_1/stack�
Fsequential_4/sequential_3/random_translation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2H
Fsequential_4/sequential_3/random_translation_1/strided_slice_1/stack_1�
Fsequential_4/sequential_3/random_translation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2H
Fsequential_4/sequential_3/random_translation_1/strided_slice_1/stack_2�
>sequential_4/sequential_3/random_translation_1/strided_slice_1StridedSlice=sequential_4/sequential_3/random_translation_1/Shape:output:0Msequential_4/sequential_3/random_translation_1/strided_slice_1/stack:output:0Osequential_4/sequential_3/random_translation_1/strided_slice_1/stack_1:output:0Osequential_4/sequential_3/random_translation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2@
>sequential_4/sequential_3/random_translation_1/strided_slice_1�
3sequential_4/sequential_3/random_translation_1/CastCastGsequential_4/sequential_3/random_translation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 25
3sequential_4/sequential_3/random_translation_1/Cast�
Dsequential_4/sequential_3/random_translation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2F
Dsequential_4/sequential_3/random_translation_1/strided_slice_2/stack�
Fsequential_4/sequential_3/random_translation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2H
Fsequential_4/sequential_3/random_translation_1/strided_slice_2/stack_1�
Fsequential_4/sequential_3/random_translation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2H
Fsequential_4/sequential_3/random_translation_1/strided_slice_2/stack_2�
>sequential_4/sequential_3/random_translation_1/strided_slice_2StridedSlice=sequential_4/sequential_3/random_translation_1/Shape:output:0Msequential_4/sequential_3/random_translation_1/strided_slice_2/stack:output:0Osequential_4/sequential_3/random_translation_1/strided_slice_2/stack_1:output:0Osequential_4/sequential_3/random_translation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2@
>sequential_4/sequential_3/random_translation_1/strided_slice_2�
5sequential_4/sequential_3/random_translation_1/Cast_1CastGsequential_4/sequential_3/random_translation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 27
5sequential_4/sequential_3/random_translation_1/Cast_1�
Gsequential_4/sequential_3/random_translation_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2I
Gsequential_4/sequential_3/random_translation_1/stateful_uniform/shape/1�
Esequential_4/sequential_3/random_translation_1/stateful_uniform/shapePackEsequential_4/sequential_3/random_translation_1/strided_slice:output:0Psequential_4/sequential_3/random_translation_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2G
Esequential_4/sequential_3/random_translation_1/stateful_uniform/shape�
Csequential_4/sequential_3/random_translation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2E
Csequential_4/sequential_3/random_translation_1/stateful_uniform/min�
Csequential_4/sequential_3/random_translation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *��>2E
Csequential_4/sequential_3/random_translation_1/stateful_uniform/max�
Ysequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2[
Ysequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniform/algorithm�
Osequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniformStatefulUniformXsequential_4_sequential_3_random_translation_1_stateful_uniform_statefuluniform_resourcebsequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniform/algorithm:output:0Nsequential_4/sequential_3/random_translation_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype02Q
Osequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniform�
Csequential_4/sequential_3/random_translation_1/stateful_uniform/subSubLsequential_4/sequential_3/random_translation_1/stateful_uniform/max:output:0Lsequential_4/sequential_3/random_translation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2E
Csequential_4/sequential_3/random_translation_1/stateful_uniform/sub�
Csequential_4/sequential_3/random_translation_1/stateful_uniform/mulMulXsequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniform:output:0Gsequential_4/sequential_3/random_translation_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������2E
Csequential_4/sequential_3/random_translation_1/stateful_uniform/mul�
?sequential_4/sequential_3/random_translation_1/stateful_uniformAddGsequential_4/sequential_3/random_translation_1/stateful_uniform/mul:z:0Lsequential_4/sequential_3/random_translation_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2A
?sequential_4/sequential_3/random_translation_1/stateful_uniform�
2sequential_4/sequential_3/random_translation_1/mulMulCsequential_4/sequential_3/random_translation_1/stateful_uniform:z:07sequential_4/sequential_3/random_translation_1/Cast:y:0*
T0*'
_output_shapes
:���������24
2sequential_4/sequential_3/random_translation_1/mul�
Isequential_4/sequential_3/random_translation_1/stateful_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2K
Isequential_4/sequential_3/random_translation_1/stateful_uniform_1/shape/1�
Gsequential_4/sequential_3/random_translation_1/stateful_uniform_1/shapePackEsequential_4/sequential_3/random_translation_1/strided_slice:output:0Rsequential_4/sequential_3/random_translation_1/stateful_uniform_1/shape/1:output:0*
N*
T0*
_output_shapes
:2I
Gsequential_4/sequential_3/random_translation_1/stateful_uniform_1/shape�
Esequential_4/sequential_3/random_translation_1/stateful_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2G
Esequential_4/sequential_3/random_translation_1/stateful_uniform_1/min�
Esequential_4/sequential_3/random_translation_1/stateful_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2G
Esequential_4/sequential_3/random_translation_1/stateful_uniform_1/max�
[sequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2]
[sequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniform/algorithm�
Qsequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniformStatefulUniformXsequential_4_sequential_3_random_translation_1_stateful_uniform_statefuluniform_resourcedsequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniform/algorithm:output:0Psequential_4/sequential_3/random_translation_1/stateful_uniform_1/shape:output:0P^sequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniform*'
_output_shapes
:���������*
shape_dtype02S
Qsequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniform�
Esequential_4/sequential_3/random_translation_1/stateful_uniform_1/subSubNsequential_4/sequential_3/random_translation_1/stateful_uniform_1/max:output:0Nsequential_4/sequential_3/random_translation_1/stateful_uniform_1/min:output:0*
T0*
_output_shapes
: 2G
Esequential_4/sequential_3/random_translation_1/stateful_uniform_1/sub�
Esequential_4/sequential_3/random_translation_1/stateful_uniform_1/mulMulZsequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniform:output:0Isequential_4/sequential_3/random_translation_1/stateful_uniform_1/sub:z:0*
T0*'
_output_shapes
:���������2G
Esequential_4/sequential_3/random_translation_1/stateful_uniform_1/mul�
Asequential_4/sequential_3/random_translation_1/stateful_uniform_1AddIsequential_4/sequential_3/random_translation_1/stateful_uniform_1/mul:z:0Nsequential_4/sequential_3/random_translation_1/stateful_uniform_1/min:output:0*
T0*'
_output_shapes
:���������2C
Asequential_4/sequential_3/random_translation_1/stateful_uniform_1�
4sequential_4/sequential_3/random_translation_1/mul_1MulEsequential_4/sequential_3/random_translation_1/stateful_uniform_1:z:09sequential_4/sequential_3/random_translation_1/Cast_1:y:0*
T0*'
_output_shapes
:���������26
4sequential_4/sequential_3/random_translation_1/mul_1�
:sequential_4/sequential_3/random_translation_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2<
:sequential_4/sequential_3/random_translation_1/concat/axis�
5sequential_4/sequential_3/random_translation_1/concatConcatV28sequential_4/sequential_3/random_translation_1/mul_1:z:06sequential_4/sequential_3/random_translation_1/mul:z:0Csequential_4/sequential_3/random_translation_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������27
5sequential_4/sequential_3/random_translation_1/concat�
Gsequential_4/sequential_3/random_translation_1/translation_matrix/ShapeShape>sequential_4/sequential_3/random_translation_1/concat:output:0*
T0*
_output_shapes
:2I
Gsequential_4/sequential_3/random_translation_1/translation_matrix/Shape�
Usequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2W
Usequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stack�
Wsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2Y
Wsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stack_1�
Wsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2Y
Wsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stack_2�
Osequential_4/sequential_3/random_translation_1/translation_matrix/strided_sliceStridedSlicePsequential_4/sequential_3/random_translation_1/translation_matrix/Shape:output:0^sequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stack:output:0`sequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stack_1:output:0`sequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2Q
Osequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice�
Lsequential_4/sequential_3/random_translation_1/translation_matrix/ones/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2N
Lsequential_4/sequential_3/random_translation_1/translation_matrix/ones/mul/y�
Jsequential_4/sequential_3/random_translation_1/translation_matrix/ones/mulMulXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0Usequential_4/sequential_3/random_translation_1/translation_matrix/ones/mul/y:output:0*
T0*
_output_shapes
: 2L
Jsequential_4/sequential_3/random_translation_1/translation_matrix/ones/mul�
Msequential_4/sequential_3/random_translation_1/translation_matrix/ones/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2O
Msequential_4/sequential_3/random_translation_1/translation_matrix/ones/Less/y�
Ksequential_4/sequential_3/random_translation_1/translation_matrix/ones/LessLessNsequential_4/sequential_3/random_translation_1/translation_matrix/ones/mul:z:0Vsequential_4/sequential_3/random_translation_1/translation_matrix/ones/Less/y:output:0*
T0*
_output_shapes
: 2M
Ksequential_4/sequential_3/random_translation_1/translation_matrix/ones/Less�
Osequential_4/sequential_3/random_translation_1/translation_matrix/ones/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2Q
Osequential_4/sequential_3/random_translation_1/translation_matrix/ones/packed/1�
Msequential_4/sequential_3/random_translation_1/translation_matrix/ones/packedPackXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0Xsequential_4/sequential_3/random_translation_1/translation_matrix/ones/packed/1:output:0*
N*
T0*
_output_shapes
:2O
Msequential_4/sequential_3/random_translation_1/translation_matrix/ones/packed�
Lsequential_4/sequential_3/random_translation_1/translation_matrix/ones/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2N
Lsequential_4/sequential_3/random_translation_1/translation_matrix/ones/Const�
Fsequential_4/sequential_3/random_translation_1/translation_matrix/onesFillVsequential_4/sequential_3/random_translation_1/translation_matrix/ones/packed:output:0Usequential_4/sequential_3/random_translation_1/translation_matrix/ones/Const:output:0*
T0*'
_output_shapes
:���������2H
Fsequential_4/sequential_3/random_translation_1/translation_matrix/ones�
Msequential_4/sequential_3/random_translation_1/translation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2O
Msequential_4/sequential_3/random_translation_1/translation_matrix/zeros/mul/y�
Ksequential_4/sequential_3/random_translation_1/translation_matrix/zeros/mulMulXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0Vsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2M
Ksequential_4/sequential_3/random_translation_1/translation_matrix/zeros/mul�
Nsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2P
Nsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/Less/y�
Lsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/LessLessOsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/mul:z:0Wsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2N
Lsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/Less�
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2R
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros/packed/1�
Nsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/packedPackXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0Ysequential_4/sequential_3/random_translation_1/translation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2P
Nsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/packed�
Msequential_4/sequential_3/random_translation_1/translation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2O
Msequential_4/sequential_3/random_translation_1/translation_matrix/zeros/Const�
Gsequential_4/sequential_3/random_translation_1/translation_matrix/zerosFillWsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/packed:output:0Vsequential_4/sequential_3/random_translation_1/translation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2I
Gsequential_4/sequential_3/random_translation_1/translation_matrix/zeros�
Wsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2Y
Wsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stack�
Ysequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2[
Ysequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_1�
Ysequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2[
Ysequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_2�
Qsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1StridedSlice>sequential_4/sequential_3/random_translation_1/concat:output:0`sequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stack:output:0bsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_1:output:0bsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2S
Qsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1�
Esequential_4/sequential_3/random_translation_1/translation_matrix/NegNegZsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2G
Esequential_4/sequential_3/random_translation_1/translation_matrix/Neg�
Osequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2Q
Osequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/mul/y�
Msequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/mulMulXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0Xsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2O
Msequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/mul�
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2R
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/Less/y�
Nsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/LessLessQsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/mul:z:0Ysequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2P
Nsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/Less�
Rsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2T
Rsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/packed/1�
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/packedPackXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0[sequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2R
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/packed�
Osequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2Q
Osequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/Const�
Isequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1FillYsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/packed:output:0Xsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������2K
Isequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1�
Nsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2P
Nsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/mul/y�
Lsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/mulMulXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0Wsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/mul/y:output:0*
T0*
_output_shapes
: 2N
Lsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/mul�
Osequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2Q
Osequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/Less/y�
Msequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/LessLessPsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/mul:z:0Xsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/Less/y:output:0*
T0*
_output_shapes
: 2O
Msequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/Less�
Qsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2S
Qsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/packed/1�
Osequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/packedPackXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0Zsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/packed/1:output:0*
N*
T0*
_output_shapes
:2Q
Osequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/packed�
Nsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2P
Nsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/Const�
Hsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1FillXsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/packed:output:0Wsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1/Const:output:0*
T0*'
_output_shapes
:���������2J
Hsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1�
Wsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2Y
Wsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stack�
Ysequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2[
Ysequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_1�
Ysequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2[
Ysequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_2�
Qsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2StridedSlice>sequential_4/sequential_3/random_translation_1/concat:output:0`sequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stack:output:0bsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_1:output:0bsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2S
Qsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2�
Gsequential_4/sequential_3/random_translation_1/translation_matrix/Neg_1NegZsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2I
Gsequential_4/sequential_3/random_translation_1/translation_matrix/Neg_1�
Osequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2Q
Osequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/mul/y�
Msequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/mulMulXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0Xsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2O
Msequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/mul�
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2R
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/Less/y�
Nsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/LessLessQsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/mul:z:0Ysequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2P
Nsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/Less�
Rsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2T
Rsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/packed/1�
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/packedPackXsequential_4/sequential_3/random_translation_1/translation_matrix/strided_slice:output:0[sequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2R
Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/packed�
Osequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2Q
Osequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/Const�
Isequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2FillYsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/packed:output:0Xsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������2K
Isequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2�
Msequential_4/sequential_3/random_translation_1/translation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2O
Msequential_4/sequential_3/random_translation_1/translation_matrix/concat/axis�
Hsequential_4/sequential_3/random_translation_1/translation_matrix/concatConcatV2Osequential_4/sequential_3/random_translation_1/translation_matrix/ones:output:0Psequential_4/sequential_3/random_translation_1/translation_matrix/zeros:output:0Isequential_4/sequential_3/random_translation_1/translation_matrix/Neg:y:0Rsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_1:output:0Qsequential_4/sequential_3/random_translation_1/translation_matrix/ones_1:output:0Ksequential_4/sequential_3/random_translation_1/translation_matrix/Neg_1:y:0Rsequential_4/sequential_3/random_translation_1/translation_matrix/zeros_2:output:0Vsequential_4/sequential_3/random_translation_1/translation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2J
Hsequential_4/sequential_3/random_translation_1/translation_matrix/concat�
>sequential_4/sequential_3/random_translation_1/transform/ShapeShapeesequential_4/sequential_3/random_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2@
>sequential_4/sequential_3/random_translation_1/transform/Shape�
Lsequential_4/sequential_3/random_translation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2N
Lsequential_4/sequential_3/random_translation_1/transform/strided_slice/stack�
Nsequential_4/sequential_3/random_translation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2P
Nsequential_4/sequential_3/random_translation_1/transform/strided_slice/stack_1�
Nsequential_4/sequential_3/random_translation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2P
Nsequential_4/sequential_3/random_translation_1/transform/strided_slice/stack_2�
Fsequential_4/sequential_3/random_translation_1/transform/strided_sliceStridedSliceGsequential_4/sequential_3/random_translation_1/transform/Shape:output:0Usequential_4/sequential_3/random_translation_1/transform/strided_slice/stack:output:0Wsequential_4/sequential_3/random_translation_1/transform/strided_slice/stack_1:output:0Wsequential_4/sequential_3/random_translation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2H
Fsequential_4/sequential_3/random_translation_1/transform/strided_slice�
Csequential_4/sequential_3/random_translation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2E
Csequential_4/sequential_3/random_translation_1/transform/fill_value�
Ssequential_4/sequential_3/random_translation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3esequential_4/sequential_3/random_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0Qsequential_4/sequential_3/random_translation_1/translation_matrix/concat:output:0Osequential_4/sequential_3/random_translation_1/transform/strided_slice:output:0Lsequential_4/sequential_3/random_translation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	NEAREST*
interpolation
BILINEAR2U
Ssequential_4/sequential_3/random_translation_1/transform/ImageProjectiveTransformV3�
-sequential_4/sequential_3/random_zoom_1/ShapeShapehsequential_4/sequential_3/random_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2/
-sequential_4/sequential_3/random_zoom_1/Shape�
;sequential_4/sequential_3/random_zoom_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2=
;sequential_4/sequential_3/random_zoom_1/strided_slice/stack�
=sequential_4/sequential_3/random_zoom_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2?
=sequential_4/sequential_3/random_zoom_1/strided_slice/stack_1�
=sequential_4/sequential_3/random_zoom_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2?
=sequential_4/sequential_3/random_zoom_1/strided_slice/stack_2�
5sequential_4/sequential_3/random_zoom_1/strided_sliceStridedSlice6sequential_4/sequential_3/random_zoom_1/Shape:output:0Dsequential_4/sequential_3/random_zoom_1/strided_slice/stack:output:0Fsequential_4/sequential_3/random_zoom_1/strided_slice/stack_1:output:0Fsequential_4/sequential_3/random_zoom_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask27
5sequential_4/sequential_3/random_zoom_1/strided_slice�
=sequential_4/sequential_3/random_zoom_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2?
=sequential_4/sequential_3/random_zoom_1/strided_slice_1/stack�
?sequential_4/sequential_3/random_zoom_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2A
?sequential_4/sequential_3/random_zoom_1/strided_slice_1/stack_1�
?sequential_4/sequential_3/random_zoom_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2A
?sequential_4/sequential_3/random_zoom_1/strided_slice_1/stack_2�
7sequential_4/sequential_3/random_zoom_1/strided_slice_1StridedSlice6sequential_4/sequential_3/random_zoom_1/Shape:output:0Fsequential_4/sequential_3/random_zoom_1/strided_slice_1/stack:output:0Hsequential_4/sequential_3/random_zoom_1/strided_slice_1/stack_1:output:0Hsequential_4/sequential_3/random_zoom_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask29
7sequential_4/sequential_3/random_zoom_1/strided_slice_1�
,sequential_4/sequential_3/random_zoom_1/CastCast@sequential_4/sequential_3/random_zoom_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2.
,sequential_4/sequential_3/random_zoom_1/Cast�
=sequential_4/sequential_3/random_zoom_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2?
=sequential_4/sequential_3/random_zoom_1/strided_slice_2/stack�
?sequential_4/sequential_3/random_zoom_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2A
?sequential_4/sequential_3/random_zoom_1/strided_slice_2/stack_1�
?sequential_4/sequential_3/random_zoom_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2A
?sequential_4/sequential_3/random_zoom_1/strided_slice_2/stack_2�
7sequential_4/sequential_3/random_zoom_1/strided_slice_2StridedSlice6sequential_4/sequential_3/random_zoom_1/Shape:output:0Fsequential_4/sequential_3/random_zoom_1/strided_slice_2/stack:output:0Hsequential_4/sequential_3/random_zoom_1/strided_slice_2/stack_1:output:0Hsequential_4/sequential_3/random_zoom_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask29
7sequential_4/sequential_3/random_zoom_1/strided_slice_2�
.sequential_4/sequential_3/random_zoom_1/Cast_1Cast@sequential_4/sequential_3/random_zoom_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 20
.sequential_4/sequential_3/random_zoom_1/Cast_1�
@sequential_4/sequential_3/random_zoom_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2B
@sequential_4/sequential_3/random_zoom_1/stateful_uniform/shape/1�
>sequential_4/sequential_3/random_zoom_1/stateful_uniform/shapePack>sequential_4/sequential_3/random_zoom_1/strided_slice:output:0Isequential_4/sequential_3/random_zoom_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2@
>sequential_4/sequential_3/random_zoom_1/stateful_uniform/shape�
<sequential_4/sequential_3/random_zoom_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *��L?2>
<sequential_4/sequential_3/random_zoom_1/stateful_uniform/min�
<sequential_4/sequential_3/random_zoom_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���?2>
<sequential_4/sequential_3/random_zoom_1/stateful_uniform/max�
Rsequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2T
Rsequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniform/algorithm�
Hsequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniformStatefulUniformQsequential_4_sequential_3_random_zoom_1_stateful_uniform_statefuluniform_resource[sequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniform/algorithm:output:0Gsequential_4/sequential_3/random_zoom_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype02J
Hsequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniform�
<sequential_4/sequential_3/random_zoom_1/stateful_uniform/subSubEsequential_4/sequential_3/random_zoom_1/stateful_uniform/max:output:0Esequential_4/sequential_3/random_zoom_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2>
<sequential_4/sequential_3/random_zoom_1/stateful_uniform/sub�
<sequential_4/sequential_3/random_zoom_1/stateful_uniform/mulMulQsequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniform:output:0@sequential_4/sequential_3/random_zoom_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������2>
<sequential_4/sequential_3/random_zoom_1/stateful_uniform/mul�
8sequential_4/sequential_3/random_zoom_1/stateful_uniformAdd@sequential_4/sequential_3/random_zoom_1/stateful_uniform/mul:z:0Esequential_4/sequential_3/random_zoom_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2:
8sequential_4/sequential_3/random_zoom_1/stateful_uniform�
3sequential_4/sequential_3/random_zoom_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :25
3sequential_4/sequential_3/random_zoom_1/concat/axis�
.sequential_4/sequential_3/random_zoom_1/concatConcatV2<sequential_4/sequential_3/random_zoom_1/stateful_uniform:z:0<sequential_4/sequential_3/random_zoom_1/stateful_uniform:z:0<sequential_4/sequential_3/random_zoom_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������20
.sequential_4/sequential_3/random_zoom_1/concat�
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/ShapeShape7sequential_4/sequential_3/random_zoom_1/concat:output:0*
T0*
_output_shapes
:2;
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/Shape�
Gsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2I
Gsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack�
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2K
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_1�
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2K
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_2�
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_sliceStridedSliceBsequential_4/sequential_3/random_zoom_1/zoom_matrix/Shape:output:0Psequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack:output:0Rsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_1:output:0Rsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2C
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice�
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2;
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub/y�
7sequential_4/sequential_3/random_zoom_1/zoom_matrix/subSub2sequential_4/sequential_3/random_zoom_1/Cast_1:y:0Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 29
7sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub�
=sequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2?
=sequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv/y�
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/truedivRealDiv;sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub:z:0Fsequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 2=
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv�
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2K
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack�
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2M
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_1�
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2M
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_2�
Csequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1StridedSlice7sequential_4/sequential_3/random_zoom_1/concat:output:0Rsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack:output:0Tsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_1:output:0Tsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2E
Csequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1�
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2=
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_1/x�
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_1SubDsequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_1/x:output:0Lsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2;
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_1�
7sequential_4/sequential_3/random_zoom_1/zoom_matrix/mulMul?sequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv:z:0=sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:���������29
7sequential_4/sequential_3/random_zoom_1/zoom_matrix/mul�
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2=
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_2/y�
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_2Sub0sequential_4/sequential_3/random_zoom_1/Cast:y:0Dsequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2;
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_2�
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2A
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv_1/y�
=sequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv_1RealDiv=sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_2:z:0Hsequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 2?
=sequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv_1�
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2K
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack�
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2M
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_1�
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2M
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_2�
Csequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2StridedSlice7sequential_4/sequential_3/random_zoom_1/concat:output:0Rsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack:output:0Tsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_1:output:0Tsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2E
Csequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2�
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2=
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_3/x�
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_3SubDsequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_3/x:output:0Lsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2;
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_3�
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/mul_1MulAsequential_4/sequential_3/random_zoom_1/zoom_matrix/truediv_1:z:0=sequential_4/sequential_3/random_zoom_1/zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:���������2;
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/mul_1�
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2K
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack�
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2M
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_1�
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2M
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_2�
Csequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3StridedSlice7sequential_4/sequential_3/random_zoom_1/concat:output:0Rsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack:output:0Tsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_1:output:0Tsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2E
Csequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3�
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2A
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/mul/y�
=sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/mulMulJsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0Hsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2?
=sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/mul�
@sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2B
@sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/Less/y�
>sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/LessLessAsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/mul:z:0Isequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2@
>sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/Less�
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2D
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/packed/1�
@sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/packedPackJsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2B
@sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/packed�
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2A
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/Const�
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/zerosFillIsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/packed:output:0Hsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2;
9sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros�
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2C
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul/y�
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/mulMulJsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0Jsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2A
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul�
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2D
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/Less/y�
@sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/LessLessCsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul:z:0Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2B
@sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/Less�
Dsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2F
Dsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed/1�
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/packedPackJsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0Msequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2D
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed�
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2C
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/Const�
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1FillKsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed:output:0Jsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������2=
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1�
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2K
Isequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack�
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2M
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_1�
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2M
Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_2�
Csequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4StridedSlice7sequential_4/sequential_3/random_zoom_1/concat:output:0Rsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack:output:0Tsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_1:output:0Tsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2E
Csequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4�
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2C
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul/y�
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/mulMulJsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0Jsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2A
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul�
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2D
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/Less/y�
@sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/LessLessCsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul:z:0Ksequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2B
@sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/Less�
Dsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2F
Dsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed/1�
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/packedPackJsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0Msequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2D
Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed�
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2C
Asequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/Const�
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2FillKsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed:output:0Jsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������2=
;sequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2�
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2A
?sequential_4/sequential_3/random_zoom_1/zoom_matrix/concat/axis�
:sequential_4/sequential_3/random_zoom_1/zoom_matrix/concatConcatV2Lsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_3:output:0Bsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros:output:0;sequential_4/sequential_3/random_zoom_1/zoom_matrix/mul:z:0Dsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_1:output:0Lsequential_4/sequential_3/random_zoom_1/zoom_matrix/strided_slice_4:output:0=sequential_4/sequential_3/random_zoom_1/zoom_matrix/mul_1:z:0Dsequential_4/sequential_3/random_zoom_1/zoom_matrix/zeros_2:output:0Hsequential_4/sequential_3/random_zoom_1/zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2<
:sequential_4/sequential_3/random_zoom_1/zoom_matrix/concat�
7sequential_4/sequential_3/random_zoom_1/transform/ShapeShapehsequential_4/sequential_3/random_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:29
7sequential_4/sequential_3/random_zoom_1/transform/Shape�
Esequential_4/sequential_3/random_zoom_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2G
Esequential_4/sequential_3/random_zoom_1/transform/strided_slice/stack�
Gsequential_4/sequential_3/random_zoom_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2I
Gsequential_4/sequential_3/random_zoom_1/transform/strided_slice/stack_1�
Gsequential_4/sequential_3/random_zoom_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2I
Gsequential_4/sequential_3/random_zoom_1/transform/strided_slice/stack_2�
?sequential_4/sequential_3/random_zoom_1/transform/strided_sliceStridedSlice@sequential_4/sequential_3/random_zoom_1/transform/Shape:output:0Nsequential_4/sequential_3/random_zoom_1/transform/strided_slice/stack:output:0Psequential_4/sequential_3/random_zoom_1/transform/strided_slice/stack_1:output:0Psequential_4/sequential_3/random_zoom_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2A
?sequential_4/sequential_3/random_zoom_1/transform/strided_slice�
<sequential_4/sequential_3/random_zoom_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2>
<sequential_4/sequential_3/random_zoom_1/transform/fill_value�
Lsequential_4/sequential_3/random_zoom_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3hsequential_4/sequential_3/random_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0Csequential_4/sequential_3/random_zoom_1/zoom_matrix/concat:output:0Hsequential_4/sequential_3/random_zoom_1/transform/strided_slice:output:0Esequential_4/sequential_3/random_zoom_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2N
Lsequential_4/sequential_3/random_zoom_1/transform/ImageProjectiveTransformV3�
sequential_4/rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2!
sequential_4/rescaling_3/Cast/x�
!sequential_4/rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!sequential_4/rescaling_3/Cast_1/x�
sequential_4/rescaling_3/mulMulasequential_4/sequential_3/random_zoom_1/transform/ImageProjectiveTransformV3:transformed_images:0(sequential_4/rescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2
sequential_4/rescaling_3/mul�
sequential_4/rescaling_3/addAddV2 sequential_4/rescaling_3/mul:z:0*sequential_4/rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2
sequential_4/rescaling_3/add�
+sequential_4/conv2d_3/Conv2D/ReadVariableOpReadVariableOp4sequential_4_conv2d_3_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02-
+sequential_4/conv2d_3/Conv2D/ReadVariableOp�
sequential_4/conv2d_3/Conv2DConv2D sequential_4/rescaling_3/add:z:03sequential_4/conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&*
paddingSAME*
strides
2
sequential_4/conv2d_3/Conv2D�
,sequential_4/conv2d_3/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_conv2d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,sequential_4/conv2d_3/BiasAdd/ReadVariableOp�
sequential_4/conv2d_3/BiasAddBiasAdd%sequential_4/conv2d_3/Conv2D:output:04sequential_4/conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&2
sequential_4/conv2d_3/BiasAdd�
sequential_4/conv2d_3/ReluRelu&sequential_4/conv2d_3/BiasAdd:output:0*
T0*/
_output_shapes
:���������&2
sequential_4/conv2d_3/Relu�
$sequential_4/dropout_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2&
$sequential_4/dropout_4/dropout/Const�
"sequential_4/dropout_4/dropout/MulMul(sequential_4/conv2d_3/Relu:activations:0-sequential_4/dropout_4/dropout/Const:output:0*
T0*/
_output_shapes
:���������&2$
"sequential_4/dropout_4/dropout/Mul�
$sequential_4/dropout_4/dropout/ShapeShape(sequential_4/conv2d_3/Relu:activations:0*
T0*
_output_shapes
:2&
$sequential_4/dropout_4/dropout/Shape�
;sequential_4/dropout_4/dropout/random_uniform/RandomUniformRandomUniform-sequential_4/dropout_4/dropout/Shape:output:0*
T0*/
_output_shapes
:���������&*
dtype02=
;sequential_4/dropout_4/dropout/random_uniform/RandomUniform�
-sequential_4/dropout_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2/
-sequential_4/dropout_4/dropout/GreaterEqual/y�
+sequential_4/dropout_4/dropout/GreaterEqualGreaterEqualDsequential_4/dropout_4/dropout/random_uniform/RandomUniform:output:06sequential_4/dropout_4/dropout/GreaterEqual/y:output:0*
T0*/
_output_shapes
:���������&2-
+sequential_4/dropout_4/dropout/GreaterEqual�
#sequential_4/dropout_4/dropout/CastCast/sequential_4/dropout_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*/
_output_shapes
:���������&2%
#sequential_4/dropout_4/dropout/Cast�
$sequential_4/dropout_4/dropout/Mul_1Mul&sequential_4/dropout_4/dropout/Mul:z:0'sequential_4/dropout_4/dropout/Cast:y:0*
T0*/
_output_shapes
:���������&2&
$sequential_4/dropout_4/dropout/Mul_1�
$sequential_4/max_pooling2d_3/MaxPoolMaxPool(sequential_4/dropout_4/dropout/Mul_1:z:0*/
_output_shapes
:���������*
ksize
*
paddingVALID*
strides
2&
$sequential_4/max_pooling2d_3/MaxPool�
+sequential_4/conv2d_4/Conv2D/ReadVariableOpReadVariableOp4sequential_4_conv2d_4_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02-
+sequential_4/conv2d_4/Conv2D/ReadVariableOp�
sequential_4/conv2d_4/Conv2DConv2D-sequential_4/max_pooling2d_3/MaxPool:output:03sequential_4/conv2d_4/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������*
paddingSAME*
strides
2
sequential_4/conv2d_4/Conv2D�
,sequential_4/conv2d_4/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_conv2d_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,sequential_4/conv2d_4/BiasAdd/ReadVariableOp�
sequential_4/conv2d_4/BiasAddBiasAdd%sequential_4/conv2d_4/Conv2D:output:04sequential_4/conv2d_4/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������2
sequential_4/conv2d_4/BiasAdd�
sequential_4/conv2d_4/ReluRelu&sequential_4/conv2d_4/BiasAdd:output:0*
T0*/
_output_shapes
:���������2
sequential_4/conv2d_4/Relu�
$sequential_4/dropout_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2&
$sequential_4/dropout_5/dropout/Const�
"sequential_4/dropout_5/dropout/MulMul(sequential_4/conv2d_4/Relu:activations:0-sequential_4/dropout_5/dropout/Const:output:0*
T0*/
_output_shapes
:���������2$
"sequential_4/dropout_5/dropout/Mul�
$sequential_4/dropout_5/dropout/ShapeShape(sequential_4/conv2d_4/Relu:activations:0*
T0*
_output_shapes
:2&
$sequential_4/dropout_5/dropout/Shape�
;sequential_4/dropout_5/dropout/random_uniform/RandomUniformRandomUniform-sequential_4/dropout_5/dropout/Shape:output:0*
T0*/
_output_shapes
:���������*
dtype02=
;sequential_4/dropout_5/dropout/random_uniform/RandomUniform�
-sequential_4/dropout_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2/
-sequential_4/dropout_5/dropout/GreaterEqual/y�
+sequential_4/dropout_5/dropout/GreaterEqualGreaterEqualDsequential_4/dropout_5/dropout/random_uniform/RandomUniform:output:06sequential_4/dropout_5/dropout/GreaterEqual/y:output:0*
T0*/
_output_shapes
:���������2-
+sequential_4/dropout_5/dropout/GreaterEqual�
#sequential_4/dropout_5/dropout/CastCast/sequential_4/dropout_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*/
_output_shapes
:���������2%
#sequential_4/dropout_5/dropout/Cast�
$sequential_4/dropout_5/dropout/Mul_1Mul&sequential_4/dropout_5/dropout/Mul:z:0'sequential_4/dropout_5/dropout/Cast:y:0*
T0*/
_output_shapes
:���������2&
$sequential_4/dropout_5/dropout/Mul_1�
$sequential_4/max_pooling2d_4/MaxPoolMaxPool(sequential_4/dropout_5/dropout/Mul_1:z:0*/
_output_shapes
:���������	*
ksize
*
paddingVALID*
strides
2&
$sequential_4/max_pooling2d_4/MaxPool�
+sequential_4/conv2d_5/Conv2D/ReadVariableOpReadVariableOp4sequential_4_conv2d_5_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype02-
+sequential_4/conv2d_5/Conv2D/ReadVariableOp�
sequential_4/conv2d_5/Conv2DConv2D-sequential_4/max_pooling2d_4/MaxPool:output:03sequential_4/conv2d_5/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 *
paddingSAME*
strides
2
sequential_4/conv2d_5/Conv2D�
,sequential_4/conv2d_5/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_conv2d_5_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02.
,sequential_4/conv2d_5/BiasAdd/ReadVariableOp�
sequential_4/conv2d_5/BiasAddBiasAdd%sequential_4/conv2d_5/Conv2D:output:04sequential_4/conv2d_5/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 2
sequential_4/conv2d_5/BiasAdd�
sequential_4/conv2d_5/ReluRelu&sequential_4/conv2d_5/BiasAdd:output:0*
T0*/
_output_shapes
:���������	 2
sequential_4/conv2d_5/Relu�
$sequential_4/max_pooling2d_5/MaxPoolMaxPool(sequential_4/conv2d_5/Relu:activations:0*/
_output_shapes
:��������� *
ksize
*
paddingVALID*
strides
2&
$sequential_4/max_pooling2d_5/MaxPool�
sequential_4/flatten_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"�����  2
sequential_4/flatten_1/Const�
sequential_4/flatten_1/ReshapeReshape-sequential_4/max_pooling2d_5/MaxPool:output:0%sequential_4/flatten_1/Const:output:0*
T0*(
_output_shapes
:����������2 
sequential_4/flatten_1/Reshape�
*sequential_4/dense_3/MatMul/ReadVariableOpReadVariableOp3sequential_4_dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype02,
*sequential_4/dense_3/MatMul/ReadVariableOp�
sequential_4/dense_3/MatMulMatMul'sequential_4/flatten_1/Reshape:output:02sequential_4/dense_3/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
sequential_4/dense_3/MatMul�
+sequential_4/dense_3/BiasAdd/ReadVariableOpReadVariableOp4sequential_4_dense_3_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype02-
+sequential_4/dense_3/BiasAdd/ReadVariableOp�
sequential_4/dense_3/BiasAddBiasAdd%sequential_4/dense_3/MatMul:product:03sequential_4/dense_3/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
sequential_4/dense_3/BiasAdd�
sequential_4/dense_3/ReluRelu%sequential_4/dense_3/BiasAdd:output:0*
T0*(
_output_shapes
:����������2
sequential_4/dense_3/Relu�
$sequential_4/dropout_6/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2&
$sequential_4/dropout_6/dropout/Const�
"sequential_4/dropout_6/dropout/MulMul'sequential_4/dense_3/Relu:activations:0-sequential_4/dropout_6/dropout/Const:output:0*
T0*(
_output_shapes
:����������2$
"sequential_4/dropout_6/dropout/Mul�
$sequential_4/dropout_6/dropout/ShapeShape'sequential_4/dense_3/Relu:activations:0*
T0*
_output_shapes
:2&
$sequential_4/dropout_6/dropout/Shape�
;sequential_4/dropout_6/dropout/random_uniform/RandomUniformRandomUniform-sequential_4/dropout_6/dropout/Shape:output:0*
T0*(
_output_shapes
:����������*
dtype02=
;sequential_4/dropout_6/dropout/random_uniform/RandomUniform�
-sequential_4/dropout_6/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2/
-sequential_4/dropout_6/dropout/GreaterEqual/y�
+sequential_4/dropout_6/dropout/GreaterEqualGreaterEqualDsequential_4/dropout_6/dropout/random_uniform/RandomUniform:output:06sequential_4/dropout_6/dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:����������2-
+sequential_4/dropout_6/dropout/GreaterEqual�
#sequential_4/dropout_6/dropout/CastCast/sequential_4/dropout_6/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:����������2%
#sequential_4/dropout_6/dropout/Cast�
$sequential_4/dropout_6/dropout/Mul_1Mul&sequential_4/dropout_6/dropout/Mul:z:0'sequential_4/dropout_6/dropout/Cast:y:0*
T0*(
_output_shapes
:����������2&
$sequential_4/dropout_6/dropout/Mul_1�
*sequential_4/dense_4/MatMul/ReadVariableOpReadVariableOp3sequential_4_dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype02,
*sequential_4/dense_4/MatMul/ReadVariableOp�
sequential_4/dense_4/MatMulMatMul(sequential_4/dropout_6/dropout/Mul_1:z:02sequential_4/dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_4/MatMul�
+sequential_4/dense_4/BiasAdd/ReadVariableOpReadVariableOp4sequential_4_dense_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+sequential_4/dense_4/BiasAdd/ReadVariableOp�
sequential_4/dense_4/BiasAddBiasAdd%sequential_4/dense_4/MatMul:product:03sequential_4/dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_4/BiasAdd�
sequential_4/dense_4/ReluRelu%sequential_4/dense_4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_4/Relu�
$sequential_4/dropout_7/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2&
$sequential_4/dropout_7/dropout/Const�
"sequential_4/dropout_7/dropout/MulMul'sequential_4/dense_4/Relu:activations:0-sequential_4/dropout_7/dropout/Const:output:0*
T0*'
_output_shapes
:���������2$
"sequential_4/dropout_7/dropout/Mul�
$sequential_4/dropout_7/dropout/ShapeShape'sequential_4/dense_4/Relu:activations:0*
T0*
_output_shapes
:2&
$sequential_4/dropout_7/dropout/Shape�
;sequential_4/dropout_7/dropout/random_uniform/RandomUniformRandomUniform-sequential_4/dropout_7/dropout/Shape:output:0*
T0*'
_output_shapes
:���������*
dtype02=
;sequential_4/dropout_7/dropout/random_uniform/RandomUniform�
-sequential_4/dropout_7/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���>2/
-sequential_4/dropout_7/dropout/GreaterEqual/y�
+sequential_4/dropout_7/dropout/GreaterEqualGreaterEqualDsequential_4/dropout_7/dropout/random_uniform/RandomUniform:output:06sequential_4/dropout_7/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:���������2-
+sequential_4/dropout_7/dropout/GreaterEqual�
#sequential_4/dropout_7/dropout/CastCast/sequential_4/dropout_7/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:���������2%
#sequential_4/dropout_7/dropout/Cast�
$sequential_4/dropout_7/dropout/Mul_1Mul&sequential_4/dropout_7/dropout/Mul:z:0'sequential_4/dropout_7/dropout/Cast:y:0*
T0*'
_output_shapes
:���������2&
$sequential_4/dropout_7/dropout/Mul_1�
*sequential_4/dense_5/MatMul/ReadVariableOpReadVariableOp3sequential_4_dense_5_matmul_readvariableop_resource*
_output_shapes

:*
dtype02,
*sequential_4/dense_5/MatMul/ReadVariableOp�
sequential_4/dense_5/MatMulMatMul(sequential_4/dropout_7/dropout/Mul_1:z:02sequential_4/dense_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_5/MatMul�
+sequential_4/dense_5/BiasAdd/ReadVariableOpReadVariableOp4sequential_4_dense_5_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+sequential_4/dense_5/BiasAdd/ReadVariableOp�
sequential_4/dense_5/BiasAddBiasAdd%sequential_4/dense_5/MatMul:product:03sequential_4/dense_5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_5/BiasAdd�
sequential_4/dense_5/SigmoidSigmoid%sequential_4/dense_5/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_5/Sigmoid�
softmax_1/SoftmaxSoftmax sequential_4/dense_5/Sigmoid:y:0*
T0*'
_output_shapes
:���������2
softmax_1/Softmax�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp3sequential_4_dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOp3sequential_4_dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentitysoftmax_1/Softmax:softmax:01^dense_3/kernel/Regularizer/Square/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp-^sequential_4/conv2d_3/BiasAdd/ReadVariableOp,^sequential_4/conv2d_3/Conv2D/ReadVariableOp-^sequential_4/conv2d_4/BiasAdd/ReadVariableOp,^sequential_4/conv2d_4/Conv2D/ReadVariableOp-^sequential_4/conv2d_5/BiasAdd/ReadVariableOp,^sequential_4/conv2d_5/Conv2D/ReadVariableOp,^sequential_4/dense_3/BiasAdd/ReadVariableOp+^sequential_4/dense_3/MatMul/ReadVariableOp,^sequential_4/dense_4/BiasAdd/ReadVariableOp+^sequential_4/dense_4/MatMul/ReadVariableOp,^sequential_4/dense_5/BiasAdd/ReadVariableOp+^sequential_4/dense_5/MatMul/ReadVariableOpM^sequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniformP^sequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniformR^sequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniformI^sequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniform*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2\
,sequential_4/conv2d_3/BiasAdd/ReadVariableOp,sequential_4/conv2d_3/BiasAdd/ReadVariableOp2Z
+sequential_4/conv2d_3/Conv2D/ReadVariableOp+sequential_4/conv2d_3/Conv2D/ReadVariableOp2\
,sequential_4/conv2d_4/BiasAdd/ReadVariableOp,sequential_4/conv2d_4/BiasAdd/ReadVariableOp2Z
+sequential_4/conv2d_4/Conv2D/ReadVariableOp+sequential_4/conv2d_4/Conv2D/ReadVariableOp2\
,sequential_4/conv2d_5/BiasAdd/ReadVariableOp,sequential_4/conv2d_5/BiasAdd/ReadVariableOp2Z
+sequential_4/conv2d_5/Conv2D/ReadVariableOp+sequential_4/conv2d_5/Conv2D/ReadVariableOp2Z
+sequential_4/dense_3/BiasAdd/ReadVariableOp+sequential_4/dense_3/BiasAdd/ReadVariableOp2X
*sequential_4/dense_3/MatMul/ReadVariableOp*sequential_4/dense_3/MatMul/ReadVariableOp2Z
+sequential_4/dense_4/BiasAdd/ReadVariableOp+sequential_4/dense_4/BiasAdd/ReadVariableOp2X
*sequential_4/dense_4/MatMul/ReadVariableOp*sequential_4/dense_4/MatMul/ReadVariableOp2Z
+sequential_4/dense_5/BiasAdd/ReadVariableOp+sequential_4/dense_5/BiasAdd/ReadVariableOp2X
*sequential_4/dense_5/MatMul/ReadVariableOp*sequential_4/dense_5/MatMul/ReadVariableOp2�
Lsequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniformLsequential_4/sequential_3/random_rotation_1/stateful_uniform/StatefulUniform2�
Osequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniformOsequential_4/sequential_3/random_translation_1/stateful_uniform/StatefulUniform2�
Qsequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniformQsequential_4/sequential_3/random_translation_1/stateful_uniform_1/StatefulUniform2�
Hsequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniformHsequential_4/sequential_3/random_zoom_1/stateful_uniform/StatefulUniform:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
�
B__inference_dense_3_layer_call_and_return_conditional_losses_33959

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�0dense_3/kernel/Regularizer/Square/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:����������2
Relu�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�"
�
G__inference_sequential_5_layer_call_and_return_conditional_losses_32117
sequential_4_input
sequential_4_32078
sequential_4_32080
sequential_4_32082
sequential_4_32084
sequential_4_32086
sequential_4_32088
sequential_4_32090
sequential_4_32092
sequential_4_32094
sequential_4_32096
sequential_4_32098
sequential_4_32100
identity��0dense_3/kernel/Regularizer/Square/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�$sequential_4/StatefulPartitionedCall�
$sequential_4/StatefulPartitionedCallStatefulPartitionedCallsequential_4_inputsequential_4_32078sequential_4_32080sequential_4_32082sequential_4_32084sequential_4_32086sequential_4_32088sequential_4_32090sequential_4_32092sequential_4_32094sequential_4_32096sequential_4_32098sequential_4_32100*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_4_layer_call_and_return_conditional_losses_319232&
$sequential_4/StatefulPartitionedCall�
softmax_1/PartitionedCallPartitionedCall-sequential_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_softmax_1_layer_call_and_return_conditional_losses_320542
softmax_1/PartitionedCall�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_32090* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_32094*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity"softmax_1/PartitionedCall:output:01^dense_3/kernel/Regularizer/Square/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp%^sequential_4/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_4_input
�
f
J__inference_max_pooling2d_4_layer_call_and_return_conditional_losses_31329

inputs
identity�
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4������������������������������������*
ksize
*
paddingVALID*
strides
2	
MaxPool�
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4������������������������������������2

Identity"
identityIdentity:output:0*I
_input_shapes8
6:4������������������������������������:r n
J
_output_shapes8
6:4������������������������������������
 
_user_specified_nameinputs
��
�
G__inference_sequential_3_layer_call_and_return_conditional_losses_33791

inputs?
;random_rotation_1_stateful_uniform_statefuluniform_resourceB
>random_translation_1_stateful_uniform_statefuluniform_resource;
7random_zoom_1_stateful_uniform_statefuluniform_resource
identity��2random_rotation_1/stateful_uniform/StatefulUniform�5random_translation_1/stateful_uniform/StatefulUniform�7random_translation_1/stateful_uniform_1/StatefulUniform�.random_zoom_1/stateful_uniform/StatefulUniform�
7random_flip_1/random_flip_left_right/control_dependencyIdentityinputs*
T0*
_class
loc:@inputs*/
_output_shapes
:���������&29
7random_flip_1/random_flip_left_right/control_dependency�
*random_flip_1/random_flip_left_right/ShapeShape@random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:2,
*random_flip_1/random_flip_left_right/Shape�
8random_flip_1/random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2:
8random_flip_1/random_flip_left_right/strided_slice/stack�
:random_flip_1/random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2<
:random_flip_1/random_flip_left_right/strided_slice/stack_1�
:random_flip_1/random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2<
:random_flip_1/random_flip_left_right/strided_slice/stack_2�
2random_flip_1/random_flip_left_right/strided_sliceStridedSlice3random_flip_1/random_flip_left_right/Shape:output:0Arandom_flip_1/random_flip_left_right/strided_slice/stack:output:0Crandom_flip_1/random_flip_left_right/strided_slice/stack_1:output:0Crandom_flip_1/random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask24
2random_flip_1/random_flip_left_right/strided_slice�
9random_flip_1/random_flip_left_right/random_uniform/shapePack;random_flip_1/random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2;
9random_flip_1/random_flip_left_right/random_uniform/shape�
7random_flip_1/random_flip_left_right/random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    29
7random_flip_1/random_flip_left_right/random_uniform/min�
7random_flip_1/random_flip_left_right/random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  �?29
7random_flip_1/random_flip_left_right/random_uniform/max�
Arandom_flip_1/random_flip_left_right/random_uniform/RandomUniformRandomUniformBrandom_flip_1/random_flip_left_right/random_uniform/shape:output:0*
T0*#
_output_shapes
:���������*
dtype02C
Arandom_flip_1/random_flip_left_right/random_uniform/RandomUniform�
7random_flip_1/random_flip_left_right/random_uniform/MulMulJrandom_flip_1/random_flip_left_right/random_uniform/RandomUniform:output:0@random_flip_1/random_flip_left_right/random_uniform/max:output:0*
T0*#
_output_shapes
:���������29
7random_flip_1/random_flip_left_right/random_uniform/Mul�
4random_flip_1/random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/1�
4random_flip_1/random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/2�
4random_flip_1/random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/3�
2random_flip_1/random_flip_left_right/Reshape/shapePack;random_flip_1/random_flip_left_right/strided_slice:output:0=random_flip_1/random_flip_left_right/Reshape/shape/1:output:0=random_flip_1/random_flip_left_right/Reshape/shape/2:output:0=random_flip_1/random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:24
2random_flip_1/random_flip_left_right/Reshape/shape�
,random_flip_1/random_flip_left_right/ReshapeReshape;random_flip_1/random_flip_left_right/random_uniform/Mul:z:0;random_flip_1/random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:���������2.
,random_flip_1/random_flip_left_right/Reshape�
*random_flip_1/random_flip_left_right/RoundRound5random_flip_1/random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:���������2,
*random_flip_1/random_flip_left_right/Round�
3random_flip_1/random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:25
3random_flip_1/random_flip_left_right/ReverseV2/axis�
.random_flip_1/random_flip_left_right/ReverseV2	ReverseV2@random_flip_1/random_flip_left_right/control_dependency:output:0<random_flip_1/random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:���������&20
.random_flip_1/random_flip_left_right/ReverseV2�
(random_flip_1/random_flip_left_right/mulMul.random_flip_1/random_flip_left_right/Round:y:07random_flip_1/random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:���������&2*
(random_flip_1/random_flip_left_right/mul�
*random_flip_1/random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2,
*random_flip_1/random_flip_left_right/sub/x�
(random_flip_1/random_flip_left_right/subSub3random_flip_1/random_flip_left_right/sub/x:output:0.random_flip_1/random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:���������2*
(random_flip_1/random_flip_left_right/sub�
*random_flip_1/random_flip_left_right/mul_1Mul,random_flip_1/random_flip_left_right/sub:z:0@random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:���������&2,
*random_flip_1/random_flip_left_right/mul_1�
(random_flip_1/random_flip_left_right/addAddV2,random_flip_1/random_flip_left_right/mul:z:0.random_flip_1/random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:���������&2*
(random_flip_1/random_flip_left_right/add�
random_rotation_1/ShapeShape,random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:2
random_rotation_1/Shape�
%random_rotation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2'
%random_rotation_1/strided_slice/stack�
'random_rotation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice/stack_1�
'random_rotation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice/stack_2�
random_rotation_1/strided_sliceStridedSlice random_rotation_1/Shape:output:0.random_rotation_1/strided_slice/stack:output:00random_rotation_1/strided_slice/stack_1:output:00random_rotation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2!
random_rotation_1/strided_slice�
'random_rotation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice_1/stack�
)random_rotation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_1/stack_1�
)random_rotation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_1/stack_2�
!random_rotation_1/strided_slice_1StridedSlice random_rotation_1/Shape:output:00random_rotation_1/strided_slice_1/stack:output:02random_rotation_1/strided_slice_1/stack_1:output:02random_rotation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2#
!random_rotation_1/strided_slice_1�
random_rotation_1/CastCast*random_rotation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_rotation_1/Cast�
'random_rotation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice_2/stack�
)random_rotation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_2/stack_1�
)random_rotation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_2/stack_2�
!random_rotation_1/strided_slice_2StridedSlice random_rotation_1/Shape:output:00random_rotation_1/strided_slice_2/stack:output:02random_rotation_1/strided_slice_2/stack_1:output:02random_rotation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2#
!random_rotation_1/strided_slice_2�
random_rotation_1/Cast_1Cast*random_rotation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_rotation_1/Cast_1�
(random_rotation_1/stateful_uniform/shapePack(random_rotation_1/strided_slice:output:0*
N*
T0*
_output_shapes
:2*
(random_rotation_1/stateful_uniform/shape�
&random_rotation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *   �2(
&random_rotation_1/stateful_uniform/min�
&random_rotation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2(
&random_rotation_1/stateful_uniform/max�
<random_rotation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2>
<random_rotation_1/stateful_uniform/StatefulUniform/algorithm�
2random_rotation_1/stateful_uniform/StatefulUniformStatefulUniform;random_rotation_1_stateful_uniform_statefuluniform_resourceErandom_rotation_1/stateful_uniform/StatefulUniform/algorithm:output:01random_rotation_1/stateful_uniform/shape:output:0*#
_output_shapes
:���������*
shape_dtype024
2random_rotation_1/stateful_uniform/StatefulUniform�
&random_rotation_1/stateful_uniform/subSub/random_rotation_1/stateful_uniform/max:output:0/random_rotation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2(
&random_rotation_1/stateful_uniform/sub�
&random_rotation_1/stateful_uniform/mulMul;random_rotation_1/stateful_uniform/StatefulUniform:output:0*random_rotation_1/stateful_uniform/sub:z:0*
T0*#
_output_shapes
:���������2(
&random_rotation_1/stateful_uniform/mul�
"random_rotation_1/stateful_uniformAdd*random_rotation_1/stateful_uniform/mul:z:0/random_rotation_1/stateful_uniform/min:output:0*
T0*#
_output_shapes
:���������2$
"random_rotation_1/stateful_uniform�
'random_rotation_1/rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2)
'random_rotation_1/rotation_matrix/sub/y�
%random_rotation_1/rotation_matrix/subSubrandom_rotation_1/Cast_1:y:00random_rotation_1/rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 2'
%random_rotation_1/rotation_matrix/sub�
%random_rotation_1/rotation_matrix/CosCos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Cos�
)random_rotation_1/rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_1/y�
'random_rotation_1/rotation_matrix/sub_1Subrandom_rotation_1/Cast_1:y:02random_rotation_1/rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_1�
%random_rotation_1/rotation_matrix/mulMul)random_rotation_1/rotation_matrix/Cos:y:0+random_rotation_1/rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/mul�
%random_rotation_1/rotation_matrix/SinSin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Sin�
)random_rotation_1/rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_2/y�
'random_rotation_1/rotation_matrix/sub_2Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_2�
'random_rotation_1/rotation_matrix/mul_1Mul)random_rotation_1/rotation_matrix/Sin:y:0+random_rotation_1/rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_1�
'random_rotation_1/rotation_matrix/sub_3Sub)random_rotation_1/rotation_matrix/mul:z:0+random_rotation_1/rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_3�
'random_rotation_1/rotation_matrix/sub_4Sub)random_rotation_1/rotation_matrix/sub:z:0+random_rotation_1/rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_4�
+random_rotation_1/rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2-
+random_rotation_1/rotation_matrix/truediv/y�
)random_rotation_1/rotation_matrix/truedivRealDiv+random_rotation_1/rotation_matrix/sub_4:z:04random_rotation_1/rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:���������2+
)random_rotation_1/rotation_matrix/truediv�
)random_rotation_1/rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_5/y�
'random_rotation_1/rotation_matrix/sub_5Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_5�
'random_rotation_1/rotation_matrix/Sin_1Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_1�
)random_rotation_1/rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_6/y�
'random_rotation_1/rotation_matrix/sub_6Subrandom_rotation_1/Cast_1:y:02random_rotation_1/rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_6�
'random_rotation_1/rotation_matrix/mul_2Mul+random_rotation_1/rotation_matrix/Sin_1:y:0+random_rotation_1/rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_2�
'random_rotation_1/rotation_matrix/Cos_1Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_1�
)random_rotation_1/rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_7/y�
'random_rotation_1/rotation_matrix/sub_7Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_7�
'random_rotation_1/rotation_matrix/mul_3Mul+random_rotation_1/rotation_matrix/Cos_1:y:0+random_rotation_1/rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_3�
%random_rotation_1/rotation_matrix/addAddV2+random_rotation_1/rotation_matrix/mul_2:z:0+random_rotation_1/rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/add�
'random_rotation_1/rotation_matrix/sub_8Sub+random_rotation_1/rotation_matrix/sub_5:z:0)random_rotation_1/rotation_matrix/add:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_8�
-random_rotation_1/rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2/
-random_rotation_1/rotation_matrix/truediv_1/y�
+random_rotation_1/rotation_matrix/truediv_1RealDiv+random_rotation_1/rotation_matrix/sub_8:z:06random_rotation_1/rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:���������2-
+random_rotation_1/rotation_matrix/truediv_1�
'random_rotation_1/rotation_matrix/ShapeShape&random_rotation_1/stateful_uniform:z:0*
T0*
_output_shapes
:2)
'random_rotation_1/rotation_matrix/Shape�
5random_rotation_1/rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 27
5random_rotation_1/rotation_matrix/strided_slice/stack�
7random_rotation_1/rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:29
7random_rotation_1/rotation_matrix/strided_slice/stack_1�
7random_rotation_1/rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:29
7random_rotation_1/rotation_matrix/strided_slice/stack_2�
/random_rotation_1/rotation_matrix/strided_sliceStridedSlice0random_rotation_1/rotation_matrix/Shape:output:0>random_rotation_1/rotation_matrix/strided_slice/stack:output:0@random_rotation_1/rotation_matrix/strided_slice/stack_1:output:0@random_rotation_1/rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask21
/random_rotation_1/rotation_matrix/strided_slice�
'random_rotation_1/rotation_matrix/Cos_2Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_2�
7random_rotation_1/rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_1/stack�
9random_rotation_1/rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_1/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_1/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_1StridedSlice+random_rotation_1/rotation_matrix/Cos_2:y:0@random_rotation_1/rotation_matrix/strided_slice_1/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_1/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_1�
'random_rotation_1/rotation_matrix/Sin_2Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_2�
7random_rotation_1/rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_2/stack�
9random_rotation_1/rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_2/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_2/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_2StridedSlice+random_rotation_1/rotation_matrix/Sin_2:y:0@random_rotation_1/rotation_matrix/strided_slice_2/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_2/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_2�
%random_rotation_1/rotation_matrix/NegNeg:random_rotation_1/rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Neg�
7random_rotation_1/rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_3/stack�
9random_rotation_1/rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_3/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_3/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_3StridedSlice-random_rotation_1/rotation_matrix/truediv:z:0@random_rotation_1/rotation_matrix/strided_slice_3/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_3/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_3�
'random_rotation_1/rotation_matrix/Sin_3Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_3�
7random_rotation_1/rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_4/stack�
9random_rotation_1/rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_4/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_4/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_4StridedSlice+random_rotation_1/rotation_matrix/Sin_3:y:0@random_rotation_1/rotation_matrix/strided_slice_4/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_4/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_4�
'random_rotation_1/rotation_matrix/Cos_3Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_3�
7random_rotation_1/rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_5/stack�
9random_rotation_1/rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_5/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_5/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_5StridedSlice+random_rotation_1/rotation_matrix/Cos_3:y:0@random_rotation_1/rotation_matrix/strided_slice_5/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_5/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_5�
7random_rotation_1/rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_6/stack�
9random_rotation_1/rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_6/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_6/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_6StridedSlice/random_rotation_1/rotation_matrix/truediv_1:z:0@random_rotation_1/rotation_matrix/strided_slice_6/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_6/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_6�
-random_rotation_1/rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2/
-random_rotation_1/rotation_matrix/zeros/mul/y�
+random_rotation_1/rotation_matrix/zeros/mulMul8random_rotation_1/rotation_matrix/strided_slice:output:06random_rotation_1/rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2-
+random_rotation_1/rotation_matrix/zeros/mul�
.random_rotation_1/rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�20
.random_rotation_1/rotation_matrix/zeros/Less/y�
,random_rotation_1/rotation_matrix/zeros/LessLess/random_rotation_1/rotation_matrix/zeros/mul:z:07random_rotation_1/rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2.
,random_rotation_1/rotation_matrix/zeros/Less�
0random_rotation_1/rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :22
0random_rotation_1/rotation_matrix/zeros/packed/1�
.random_rotation_1/rotation_matrix/zeros/packedPack8random_rotation_1/rotation_matrix/strided_slice:output:09random_rotation_1/rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:20
.random_rotation_1/rotation_matrix/zeros/packed�
-random_rotation_1/rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2/
-random_rotation_1/rotation_matrix/zeros/Const�
'random_rotation_1/rotation_matrix/zerosFill7random_rotation_1/rotation_matrix/zeros/packed:output:06random_rotation_1/rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/zeros�
-random_rotation_1/rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2/
-random_rotation_1/rotation_matrix/concat/axis�
(random_rotation_1/rotation_matrix/concatConcatV2:random_rotation_1/rotation_matrix/strided_slice_1:output:0)random_rotation_1/rotation_matrix/Neg:y:0:random_rotation_1/rotation_matrix/strided_slice_3:output:0:random_rotation_1/rotation_matrix/strided_slice_4:output:0:random_rotation_1/rotation_matrix/strided_slice_5:output:0:random_rotation_1/rotation_matrix/strided_slice_6:output:00random_rotation_1/rotation_matrix/zeros:output:06random_rotation_1/rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2*
(random_rotation_1/rotation_matrix/concat�
!random_rotation_1/transform/ShapeShape,random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:2#
!random_rotation_1/transform/Shape�
/random_rotation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:21
/random_rotation_1/transform/strided_slice/stack�
1random_rotation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:23
1random_rotation_1/transform/strided_slice/stack_1�
1random_rotation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:23
1random_rotation_1/transform/strided_slice/stack_2�
)random_rotation_1/transform/strided_sliceStridedSlice*random_rotation_1/transform/Shape:output:08random_rotation_1/transform/strided_slice/stack:output:0:random_rotation_1/transform/strided_slice/stack_1:output:0:random_rotation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2+
)random_rotation_1/transform/strided_slice�
&random_rotation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2(
&random_rotation_1/transform/fill_value�
6random_rotation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3,random_flip_1/random_flip_left_right/add:z:01random_rotation_1/rotation_matrix/concat:output:02random_rotation_1/transform/strided_slice:output:0/random_rotation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR28
6random_rotation_1/transform/ImageProjectiveTransformV3�
random_translation_1/ShapeShapeKrandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_translation_1/Shape�
(random_translation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2*
(random_translation_1/strided_slice/stack�
*random_translation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice/stack_1�
*random_translation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice/stack_2�
"random_translation_1/strided_sliceStridedSlice#random_translation_1/Shape:output:01random_translation_1/strided_slice/stack:output:03random_translation_1/strided_slice/stack_1:output:03random_translation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"random_translation_1/strided_slice�
*random_translation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice_1/stack�
,random_translation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_1/stack_1�
,random_translation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_1/stack_2�
$random_translation_1/strided_slice_1StridedSlice#random_translation_1/Shape:output:03random_translation_1/strided_slice_1/stack:output:05random_translation_1/strided_slice_1/stack_1:output:05random_translation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2&
$random_translation_1/strided_slice_1�
random_translation_1/CastCast-random_translation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_translation_1/Cast�
*random_translation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice_2/stack�
,random_translation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_2/stack_1�
,random_translation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_2/stack_2�
$random_translation_1/strided_slice_2StridedSlice#random_translation_1/Shape:output:03random_translation_1/strided_slice_2/stack:output:05random_translation_1/strided_slice_2/stack_1:output:05random_translation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2&
$random_translation_1/strided_slice_2�
random_translation_1/Cast_1Cast-random_translation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_translation_1/Cast_1�
-random_translation_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2/
-random_translation_1/stateful_uniform/shape/1�
+random_translation_1/stateful_uniform/shapePack+random_translation_1/strided_slice:output:06random_translation_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2-
+random_translation_1/stateful_uniform/shape�
)random_translation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2+
)random_translation_1/stateful_uniform/min�
)random_translation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *��>2+
)random_translation_1/stateful_uniform/max�
?random_translation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2A
?random_translation_1/stateful_uniform/StatefulUniform/algorithm�
5random_translation_1/stateful_uniform/StatefulUniformStatefulUniform>random_translation_1_stateful_uniform_statefuluniform_resourceHrandom_translation_1/stateful_uniform/StatefulUniform/algorithm:output:04random_translation_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype027
5random_translation_1/stateful_uniform/StatefulUniform�
)random_translation_1/stateful_uniform/subSub2random_translation_1/stateful_uniform/max:output:02random_translation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2+
)random_translation_1/stateful_uniform/sub�
)random_translation_1/stateful_uniform/mulMul>random_translation_1/stateful_uniform/StatefulUniform:output:0-random_translation_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������2+
)random_translation_1/stateful_uniform/mul�
%random_translation_1/stateful_uniformAdd-random_translation_1/stateful_uniform/mul:z:02random_translation_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2'
%random_translation_1/stateful_uniform�
random_translation_1/mulMul)random_translation_1/stateful_uniform:z:0random_translation_1/Cast:y:0*
T0*'
_output_shapes
:���������2
random_translation_1/mul�
/random_translation_1/stateful_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :21
/random_translation_1/stateful_uniform_1/shape/1�
-random_translation_1/stateful_uniform_1/shapePack+random_translation_1/strided_slice:output:08random_translation_1/stateful_uniform_1/shape/1:output:0*
N*
T0*
_output_shapes
:2/
-random_translation_1/stateful_uniform_1/shape�
+random_translation_1/stateful_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2-
+random_translation_1/stateful_uniform_1/min�
+random_translation_1/stateful_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2-
+random_translation_1/stateful_uniform_1/max�
Arandom_translation_1/stateful_uniform_1/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2C
Arandom_translation_1/stateful_uniform_1/StatefulUniform/algorithm�
7random_translation_1/stateful_uniform_1/StatefulUniformStatefulUniform>random_translation_1_stateful_uniform_statefuluniform_resourceJrandom_translation_1/stateful_uniform_1/StatefulUniform/algorithm:output:06random_translation_1/stateful_uniform_1/shape:output:06^random_translation_1/stateful_uniform/StatefulUniform*'
_output_shapes
:���������*
shape_dtype029
7random_translation_1/stateful_uniform_1/StatefulUniform�
+random_translation_1/stateful_uniform_1/subSub4random_translation_1/stateful_uniform_1/max:output:04random_translation_1/stateful_uniform_1/min:output:0*
T0*
_output_shapes
: 2-
+random_translation_1/stateful_uniform_1/sub�
+random_translation_1/stateful_uniform_1/mulMul@random_translation_1/stateful_uniform_1/StatefulUniform:output:0/random_translation_1/stateful_uniform_1/sub:z:0*
T0*'
_output_shapes
:���������2-
+random_translation_1/stateful_uniform_1/mul�
'random_translation_1/stateful_uniform_1Add/random_translation_1/stateful_uniform_1/mul:z:04random_translation_1/stateful_uniform_1/min:output:0*
T0*'
_output_shapes
:���������2)
'random_translation_1/stateful_uniform_1�
random_translation_1/mul_1Mul+random_translation_1/stateful_uniform_1:z:0random_translation_1/Cast_1:y:0*
T0*'
_output_shapes
:���������2
random_translation_1/mul_1�
 random_translation_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2"
 random_translation_1/concat/axis�
random_translation_1/concatConcatV2random_translation_1/mul_1:z:0random_translation_1/mul:z:0)random_translation_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
random_translation_1/concat�
-random_translation_1/translation_matrix/ShapeShape$random_translation_1/concat:output:0*
T0*
_output_shapes
:2/
-random_translation_1/translation_matrix/Shape�
;random_translation_1/translation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2=
;random_translation_1/translation_matrix/strided_slice/stack�
=random_translation_1/translation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2?
=random_translation_1/translation_matrix/strided_slice/stack_1�
=random_translation_1/translation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2?
=random_translation_1/translation_matrix/strided_slice/stack_2�
5random_translation_1/translation_matrix/strided_sliceStridedSlice6random_translation_1/translation_matrix/Shape:output:0Drandom_translation_1/translation_matrix/strided_slice/stack:output:0Frandom_translation_1/translation_matrix/strided_slice/stack_1:output:0Frandom_translation_1/translation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask27
5random_translation_1/translation_matrix/strided_slice�
2random_translation_1/translation_matrix/ones/mul/yConst*
_output_shapes
: *
dtype0*
value	B :24
2random_translation_1/translation_matrix/ones/mul/y�
0random_translation_1/translation_matrix/ones/mulMul>random_translation_1/translation_matrix/strided_slice:output:0;random_translation_1/translation_matrix/ones/mul/y:output:0*
T0*
_output_shapes
: 22
0random_translation_1/translation_matrix/ones/mul�
3random_translation_1/translation_matrix/ones/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�25
3random_translation_1/translation_matrix/ones/Less/y�
1random_translation_1/translation_matrix/ones/LessLess4random_translation_1/translation_matrix/ones/mul:z:0<random_translation_1/translation_matrix/ones/Less/y:output:0*
T0*
_output_shapes
: 23
1random_translation_1/translation_matrix/ones/Less�
5random_translation_1/translation_matrix/ones/packed/1Const*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/ones/packed/1�
3random_translation_1/translation_matrix/ones/packedPack>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/ones/packed/1:output:0*
N*
T0*
_output_shapes
:25
3random_translation_1/translation_matrix/ones/packed�
2random_translation_1/translation_matrix/ones/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?24
2random_translation_1/translation_matrix/ones/Const�
,random_translation_1/translation_matrix/onesFill<random_translation_1/translation_matrix/ones/packed:output:0;random_translation_1/translation_matrix/ones/Const:output:0*
T0*'
_output_shapes
:���������2.
,random_translation_1/translation_matrix/ones�
3random_translation_1/translation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :25
3random_translation_1/translation_matrix/zeros/mul/y�
1random_translation_1/translation_matrix/zeros/mulMul>random_translation_1/translation_matrix/strided_slice:output:0<random_translation_1/translation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 23
1random_translation_1/translation_matrix/zeros/mul�
4random_translation_1/translation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�26
4random_translation_1/translation_matrix/zeros/Less/y�
2random_translation_1/translation_matrix/zeros/LessLess5random_translation_1/translation_matrix/zeros/mul:z:0=random_translation_1/translation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 24
2random_translation_1/translation_matrix/zeros/Less�
6random_translation_1/translation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :28
6random_translation_1/translation_matrix/zeros/packed/1�
4random_translation_1/translation_matrix/zeros/packedPack>random_translation_1/translation_matrix/strided_slice:output:0?random_translation_1/translation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:26
4random_translation_1/translation_matrix/zeros/packed�
3random_translation_1/translation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    25
3random_translation_1/translation_matrix/zeros/Const�
-random_translation_1/translation_matrix/zerosFill=random_translation_1/translation_matrix/zeros/packed:output:0<random_translation_1/translation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2/
-random_translation_1/translation_matrix/zeros�
=random_translation_1/translation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2?
=random_translation_1/translation_matrix/strided_slice_1/stack�
?random_translation_1/translation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2A
?random_translation_1/translation_matrix/strided_slice_1/stack_1�
?random_translation_1/translation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2A
?random_translation_1/translation_matrix/strided_slice_1/stack_2�
7random_translation_1/translation_matrix/strided_slice_1StridedSlice$random_translation_1/concat:output:0Frandom_translation_1/translation_matrix/strided_slice_1/stack:output:0Hrandom_translation_1/translation_matrix/strided_slice_1/stack_1:output:0Hrandom_translation_1/translation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask29
7random_translation_1/translation_matrix/strided_slice_1�
+random_translation_1/translation_matrix/NegNeg@random_translation_1/translation_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2-
+random_translation_1/translation_matrix/Neg�
5random_translation_1/translation_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/zeros_1/mul/y�
3random_translation_1/translation_matrix/zeros_1/mulMul>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/zeros_1/mul�
6random_translation_1/translation_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�28
6random_translation_1/translation_matrix/zeros_1/Less/y�
4random_translation_1/translation_matrix/zeros_1/LessLess7random_translation_1/translation_matrix/zeros_1/mul:z:0?random_translation_1/translation_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 26
4random_translation_1/translation_matrix/zeros_1/Less�
8random_translation_1/translation_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2:
8random_translation_1/translation_matrix/zeros_1/packed/1�
6random_translation_1/translation_matrix/zeros_1/packedPack>random_translation_1/translation_matrix/strided_slice:output:0Arandom_translation_1/translation_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:28
6random_translation_1/translation_matrix/zeros_1/packed�
5random_translation_1/translation_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    27
5random_translation_1/translation_matrix/zeros_1/Const�
/random_translation_1/translation_matrix/zeros_1Fill?random_translation_1/translation_matrix/zeros_1/packed:output:0>random_translation_1/translation_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������21
/random_translation_1/translation_matrix/zeros_1�
4random_translation_1/translation_matrix/ones_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :26
4random_translation_1/translation_matrix/ones_1/mul/y�
2random_translation_1/translation_matrix/ones_1/mulMul>random_translation_1/translation_matrix/strided_slice:output:0=random_translation_1/translation_matrix/ones_1/mul/y:output:0*
T0*
_output_shapes
: 24
2random_translation_1/translation_matrix/ones_1/mul�
5random_translation_1/translation_matrix/ones_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�27
5random_translation_1/translation_matrix/ones_1/Less/y�
3random_translation_1/translation_matrix/ones_1/LessLess6random_translation_1/translation_matrix/ones_1/mul:z:0>random_translation_1/translation_matrix/ones_1/Less/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/ones_1/Less�
7random_translation_1/translation_matrix/ones_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :29
7random_translation_1/translation_matrix/ones_1/packed/1�
5random_translation_1/translation_matrix/ones_1/packedPack>random_translation_1/translation_matrix/strided_slice:output:0@random_translation_1/translation_matrix/ones_1/packed/1:output:0*
N*
T0*
_output_shapes
:27
5random_translation_1/translation_matrix/ones_1/packed�
4random_translation_1/translation_matrix/ones_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?26
4random_translation_1/translation_matrix/ones_1/Const�
.random_translation_1/translation_matrix/ones_1Fill>random_translation_1/translation_matrix/ones_1/packed:output:0=random_translation_1/translation_matrix/ones_1/Const:output:0*
T0*'
_output_shapes
:���������20
.random_translation_1/translation_matrix/ones_1�
=random_translation_1/translation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2?
=random_translation_1/translation_matrix/strided_slice_2/stack�
?random_translation_1/translation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2A
?random_translation_1/translation_matrix/strided_slice_2/stack_1�
?random_translation_1/translation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2A
?random_translation_1/translation_matrix/strided_slice_2/stack_2�
7random_translation_1/translation_matrix/strided_slice_2StridedSlice$random_translation_1/concat:output:0Frandom_translation_1/translation_matrix/strided_slice_2/stack:output:0Hrandom_translation_1/translation_matrix/strided_slice_2/stack_1:output:0Hrandom_translation_1/translation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask29
7random_translation_1/translation_matrix/strided_slice_2�
-random_translation_1/translation_matrix/Neg_1Neg@random_translation_1/translation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2/
-random_translation_1/translation_matrix/Neg_1�
5random_translation_1/translation_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/zeros_2/mul/y�
3random_translation_1/translation_matrix/zeros_2/mulMul>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/zeros_2/mul�
6random_translation_1/translation_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�28
6random_translation_1/translation_matrix/zeros_2/Less/y�
4random_translation_1/translation_matrix/zeros_2/LessLess7random_translation_1/translation_matrix/zeros_2/mul:z:0?random_translation_1/translation_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 26
4random_translation_1/translation_matrix/zeros_2/Less�
8random_translation_1/translation_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2:
8random_translation_1/translation_matrix/zeros_2/packed/1�
6random_translation_1/translation_matrix/zeros_2/packedPack>random_translation_1/translation_matrix/strided_slice:output:0Arandom_translation_1/translation_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:28
6random_translation_1/translation_matrix/zeros_2/packed�
5random_translation_1/translation_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    27
5random_translation_1/translation_matrix/zeros_2/Const�
/random_translation_1/translation_matrix/zeros_2Fill?random_translation_1/translation_matrix/zeros_2/packed:output:0>random_translation_1/translation_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������21
/random_translation_1/translation_matrix/zeros_2�
3random_translation_1/translation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :25
3random_translation_1/translation_matrix/concat/axis�
.random_translation_1/translation_matrix/concatConcatV25random_translation_1/translation_matrix/ones:output:06random_translation_1/translation_matrix/zeros:output:0/random_translation_1/translation_matrix/Neg:y:08random_translation_1/translation_matrix/zeros_1:output:07random_translation_1/translation_matrix/ones_1:output:01random_translation_1/translation_matrix/Neg_1:y:08random_translation_1/translation_matrix/zeros_2:output:0<random_translation_1/translation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������20
.random_translation_1/translation_matrix/concat�
$random_translation_1/transform/ShapeShapeKrandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2&
$random_translation_1/transform/Shape�
2random_translation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:24
2random_translation_1/transform/strided_slice/stack�
4random_translation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:26
4random_translation_1/transform/strided_slice/stack_1�
4random_translation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:26
4random_translation_1/transform/strided_slice/stack_2�
,random_translation_1/transform/strided_sliceStridedSlice-random_translation_1/transform/Shape:output:0;random_translation_1/transform/strided_slice/stack:output:0=random_translation_1/transform/strided_slice/stack_1:output:0=random_translation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2.
,random_translation_1/transform/strided_slice�
)random_translation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2+
)random_translation_1/transform/fill_value�
9random_translation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Krandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:07random_translation_1/translation_matrix/concat:output:05random_translation_1/transform/strided_slice:output:02random_translation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	NEAREST*
interpolation
BILINEAR2;
9random_translation_1/transform/ImageProjectiveTransformV3�
random_zoom_1/ShapeShapeNrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_zoom_1/Shape�
!random_zoom_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2#
!random_zoom_1/strided_slice/stack�
#random_zoom_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice/stack_1�
#random_zoom_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice/stack_2�
random_zoom_1/strided_sliceStridedSlicerandom_zoom_1/Shape:output:0*random_zoom_1/strided_slice/stack:output:0,random_zoom_1/strided_slice/stack_1:output:0,random_zoom_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice�
#random_zoom_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice_1/stack�
%random_zoom_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_1/stack_1�
%random_zoom_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_1/stack_2�
random_zoom_1/strided_slice_1StridedSlicerandom_zoom_1/Shape:output:0,random_zoom_1/strided_slice_1/stack:output:0.random_zoom_1/strided_slice_1/stack_1:output:0.random_zoom_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice_1�
random_zoom_1/CastCast&random_zoom_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_zoom_1/Cast�
#random_zoom_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice_2/stack�
%random_zoom_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_2/stack_1�
%random_zoom_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_2/stack_2�
random_zoom_1/strided_slice_2StridedSlicerandom_zoom_1/Shape:output:0,random_zoom_1/strided_slice_2/stack:output:0.random_zoom_1/strided_slice_2/stack_1:output:0.random_zoom_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice_2�
random_zoom_1/Cast_1Cast&random_zoom_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_zoom_1/Cast_1�
&random_zoom_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2(
&random_zoom_1/stateful_uniform/shape/1�
$random_zoom_1/stateful_uniform/shapePack$random_zoom_1/strided_slice:output:0/random_zoom_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2&
$random_zoom_1/stateful_uniform/shape�
"random_zoom_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *��L?2$
"random_zoom_1/stateful_uniform/min�
"random_zoom_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���?2$
"random_zoom_1/stateful_uniform/max�
8random_zoom_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2:
8random_zoom_1/stateful_uniform/StatefulUniform/algorithm�
.random_zoom_1/stateful_uniform/StatefulUniformStatefulUniform7random_zoom_1_stateful_uniform_statefuluniform_resourceArandom_zoom_1/stateful_uniform/StatefulUniform/algorithm:output:0-random_zoom_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype020
.random_zoom_1/stateful_uniform/StatefulUniform�
"random_zoom_1/stateful_uniform/subSub+random_zoom_1/stateful_uniform/max:output:0+random_zoom_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2$
"random_zoom_1/stateful_uniform/sub�
"random_zoom_1/stateful_uniform/mulMul7random_zoom_1/stateful_uniform/StatefulUniform:output:0&random_zoom_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������2$
"random_zoom_1/stateful_uniform/mul�
random_zoom_1/stateful_uniformAdd&random_zoom_1/stateful_uniform/mul:z:0+random_zoom_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2 
random_zoom_1/stateful_uniformx
random_zoom_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
random_zoom_1/concat/axis�
random_zoom_1/concatConcatV2"random_zoom_1/stateful_uniform:z:0"random_zoom_1/stateful_uniform:z:0"random_zoom_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
random_zoom_1/concat�
random_zoom_1/zoom_matrix/ShapeShaperandom_zoom_1/concat:output:0*
T0*
_output_shapes
:2!
random_zoom_1/zoom_matrix/Shape�
-random_zoom_1/zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2/
-random_zoom_1/zoom_matrix/strided_slice/stack�
/random_zoom_1/zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:21
/random_zoom_1/zoom_matrix/strided_slice/stack_1�
/random_zoom_1/zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:21
/random_zoom_1/zoom_matrix/strided_slice/stack_2�
'random_zoom_1/zoom_matrix/strided_sliceStridedSlice(random_zoom_1/zoom_matrix/Shape:output:06random_zoom_1/zoom_matrix/strided_slice/stack:output:08random_zoom_1/zoom_matrix/strided_slice/stack_1:output:08random_zoom_1/zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2)
'random_zoom_1/zoom_matrix/strided_slice�
random_zoom_1/zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2!
random_zoom_1/zoom_matrix/sub/y�
random_zoom_1/zoom_matrix/subSubrandom_zoom_1/Cast_1:y:0(random_zoom_1/zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 2
random_zoom_1/zoom_matrix/sub�
#random_zoom_1/zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2%
#random_zoom_1/zoom_matrix/truediv/y�
!random_zoom_1/zoom_matrix/truedivRealDiv!random_zoom_1/zoom_matrix/sub:z:0,random_zoom_1/zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 2#
!random_zoom_1/zoom_matrix/truediv�
/random_zoom_1/zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            21
/random_zoom_1/zoom_matrix/strided_slice_1/stack�
1random_zoom_1/zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_1/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_1/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_1StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_1/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_1/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_1�
!random_zoom_1/zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_1/x�
random_zoom_1/zoom_matrix/sub_1Sub*random_zoom_1/zoom_matrix/sub_1/x:output:02random_zoom_1/zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/sub_1�
random_zoom_1/zoom_matrix/mulMul%random_zoom_1/zoom_matrix/truediv:z:0#random_zoom_1/zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:���������2
random_zoom_1/zoom_matrix/mul�
!random_zoom_1/zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_2/y�
random_zoom_1/zoom_matrix/sub_2Subrandom_zoom_1/Cast:y:0*random_zoom_1/zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2!
random_zoom_1/zoom_matrix/sub_2�
%random_zoom_1/zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2'
%random_zoom_1/zoom_matrix/truediv_1/y�
#random_zoom_1/zoom_matrix/truediv_1RealDiv#random_zoom_1/zoom_matrix/sub_2:z:0.random_zoom_1/zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 2%
#random_zoom_1/zoom_matrix/truediv_1�
/random_zoom_1/zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom_1/zoom_matrix/strided_slice_2/stack�
1random_zoom_1/zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_2/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_2/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_2StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_2/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_2/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_2�
!random_zoom_1/zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_3/x�
random_zoom_1/zoom_matrix/sub_3Sub*random_zoom_1/zoom_matrix/sub_3/x:output:02random_zoom_1/zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/sub_3�
random_zoom_1/zoom_matrix/mul_1Mul'random_zoom_1/zoom_matrix/truediv_1:z:0#random_zoom_1/zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/mul_1�
/random_zoom_1/zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            21
/random_zoom_1/zoom_matrix/strided_slice_3/stack�
1random_zoom_1/zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_3/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_3/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_3StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_3/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_3/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_3�
%random_zoom_1/zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2'
%random_zoom_1/zoom_matrix/zeros/mul/y�
#random_zoom_1/zoom_matrix/zeros/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:0.random_zoom_1/zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2%
#random_zoom_1/zoom_matrix/zeros/mul�
&random_zoom_1/zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2(
&random_zoom_1/zoom_matrix/zeros/Less/y�
$random_zoom_1/zoom_matrix/zeros/LessLess'random_zoom_1/zoom_matrix/zeros/mul:z:0/random_zoom_1/zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2&
$random_zoom_1/zoom_matrix/zeros/Less�
(random_zoom_1/zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2*
(random_zoom_1/zoom_matrix/zeros/packed/1�
&random_zoom_1/zoom_matrix/zeros/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:01random_zoom_1/zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2(
&random_zoom_1/zoom_matrix/zeros/packed�
%random_zoom_1/zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%random_zoom_1/zoom_matrix/zeros/Const�
random_zoom_1/zoom_matrix/zerosFill/random_zoom_1/zoom_matrix/zeros/packed:output:0.random_zoom_1/zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/zeros�
'random_zoom_1/zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2)
'random_zoom_1/zoom_matrix/zeros_1/mul/y�
%random_zoom_1/zoom_matrix/zeros_1/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:00random_zoom_1/zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2'
%random_zoom_1/zoom_matrix/zeros_1/mul�
(random_zoom_1/zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2*
(random_zoom_1/zoom_matrix/zeros_1/Less/y�
&random_zoom_1/zoom_matrix/zeros_1/LessLess)random_zoom_1/zoom_matrix/zeros_1/mul:z:01random_zoom_1/zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2(
&random_zoom_1/zoom_matrix/zeros_1/Less�
*random_zoom_1/zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2,
*random_zoom_1/zoom_matrix/zeros_1/packed/1�
(random_zoom_1/zoom_matrix/zeros_1/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:03random_zoom_1/zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2*
(random_zoom_1/zoom_matrix/zeros_1/packed�
'random_zoom_1/zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2)
'random_zoom_1/zoom_matrix/zeros_1/Const�
!random_zoom_1/zoom_matrix/zeros_1Fill1random_zoom_1/zoom_matrix/zeros_1/packed:output:00random_zoom_1/zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������2#
!random_zoom_1/zoom_matrix/zeros_1�
/random_zoom_1/zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom_1/zoom_matrix/strided_slice_4/stack�
1random_zoom_1/zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_4/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_4/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_4StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_4/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_4/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_4�
'random_zoom_1/zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2)
'random_zoom_1/zoom_matrix/zeros_2/mul/y�
%random_zoom_1/zoom_matrix/zeros_2/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:00random_zoom_1/zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2'
%random_zoom_1/zoom_matrix/zeros_2/mul�
(random_zoom_1/zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2*
(random_zoom_1/zoom_matrix/zeros_2/Less/y�
&random_zoom_1/zoom_matrix/zeros_2/LessLess)random_zoom_1/zoom_matrix/zeros_2/mul:z:01random_zoom_1/zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2(
&random_zoom_1/zoom_matrix/zeros_2/Less�
*random_zoom_1/zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2,
*random_zoom_1/zoom_matrix/zeros_2/packed/1�
(random_zoom_1/zoom_matrix/zeros_2/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:03random_zoom_1/zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2*
(random_zoom_1/zoom_matrix/zeros_2/packed�
'random_zoom_1/zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2)
'random_zoom_1/zoom_matrix/zeros_2/Const�
!random_zoom_1/zoom_matrix/zeros_2Fill1random_zoom_1/zoom_matrix/zeros_2/packed:output:00random_zoom_1/zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������2#
!random_zoom_1/zoom_matrix/zeros_2�
%random_zoom_1/zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2'
%random_zoom_1/zoom_matrix/concat/axis�
 random_zoom_1/zoom_matrix/concatConcatV22random_zoom_1/zoom_matrix/strided_slice_3:output:0(random_zoom_1/zoom_matrix/zeros:output:0!random_zoom_1/zoom_matrix/mul:z:0*random_zoom_1/zoom_matrix/zeros_1:output:02random_zoom_1/zoom_matrix/strided_slice_4:output:0#random_zoom_1/zoom_matrix/mul_1:z:0*random_zoom_1/zoom_matrix/zeros_2:output:0.random_zoom_1/zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2"
 random_zoom_1/zoom_matrix/concat�
random_zoom_1/transform/ShapeShapeNrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_zoom_1/transform/Shape�
+random_zoom_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2-
+random_zoom_1/transform/strided_slice/stack�
-random_zoom_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2/
-random_zoom_1/transform/strided_slice/stack_1�
-random_zoom_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2/
-random_zoom_1/transform/strided_slice/stack_2�
%random_zoom_1/transform/strided_sliceStridedSlice&random_zoom_1/transform/Shape:output:04random_zoom_1/transform/strided_slice/stack:output:06random_zoom_1/transform/strided_slice/stack_1:output:06random_zoom_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2'
%random_zoom_1/transform/strided_slice�
"random_zoom_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2$
"random_zoom_1/transform/fill_value�
2random_zoom_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Nrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0)random_zoom_1/zoom_matrix/concat:output:0.random_zoom_1/transform/strided_slice:output:0+random_zoom_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR24
2random_zoom_1/transform/ImageProjectiveTransformV3�
IdentityIdentityGrandom_zoom_1/transform/ImageProjectiveTransformV3:transformed_images:03^random_rotation_1/stateful_uniform/StatefulUniform6^random_translation_1/stateful_uniform/StatefulUniform8^random_translation_1/stateful_uniform_1/StatefulUniform/^random_zoom_1/stateful_uniform/StatefulUniform*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*:
_input_shapes)
':���������&:::2h
2random_rotation_1/stateful_uniform/StatefulUniform2random_rotation_1/stateful_uniform/StatefulUniform2n
5random_translation_1/stateful_uniform/StatefulUniform5random_translation_1/stateful_uniform/StatefulUniform2r
7random_translation_1/stateful_uniform_1/StatefulUniform7random_translation_1/stateful_uniform_1/StatefulUniform2`
.random_zoom_1/stateful_uniform/StatefulUniform.random_zoom_1/stateful_uniform/StatefulUniform:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
c
G__inference_sequential_3_layer_call_and_return_conditional_losses_33795

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�g
�
__inference__traced_save_34286
file_prefix.
*savev2_conv2d_3_kernel_read_readvariableop,
(savev2_conv2d_3_bias_read_readvariableop.
*savev2_conv2d_4_kernel_read_readvariableop,
(savev2_conv2d_4_bias_read_readvariableop.
*savev2_conv2d_5_kernel_read_readvariableop,
(savev2_conv2d_5_bias_read_readvariableop-
)savev2_dense_3_kernel_read_readvariableop+
'savev2_dense_3_bias_read_readvariableop-
)savev2_dense_4_kernel_read_readvariableop+
'savev2_dense_4_bias_read_readvariableop-
)savev2_dense_5_kernel_read_readvariableop+
'savev2_dense_5_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop'
#savev2_variable_read_readvariableop	)
%savev2_variable_1_read_readvariableop	)
%savev2_variable_2_read_readvariableop	)
%savev2_variable_3_read_readvariableop	$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop5
1savev2_adam_conv2d_3_kernel_m_read_readvariableop3
/savev2_adam_conv2d_3_bias_m_read_readvariableop5
1savev2_adam_conv2d_4_kernel_m_read_readvariableop3
/savev2_adam_conv2d_4_bias_m_read_readvariableop5
1savev2_adam_conv2d_5_kernel_m_read_readvariableop3
/savev2_adam_conv2d_5_bias_m_read_readvariableop4
0savev2_adam_dense_3_kernel_m_read_readvariableop2
.savev2_adam_dense_3_bias_m_read_readvariableop4
0savev2_adam_dense_4_kernel_m_read_readvariableop2
.savev2_adam_dense_4_bias_m_read_readvariableop4
0savev2_adam_dense_5_kernel_m_read_readvariableop2
.savev2_adam_dense_5_bias_m_read_readvariableop5
1savev2_adam_conv2d_3_kernel_v_read_readvariableop3
/savev2_adam_conv2d_3_bias_v_read_readvariableop5
1savev2_adam_conv2d_4_kernel_v_read_readvariableop3
/savev2_adam_conv2d_4_bias_v_read_readvariableop5
1savev2_adam_conv2d_5_kernel_v_read_readvariableop3
/savev2_adam_conv2d_5_bias_v_read_readvariableop4
0savev2_adam_dense_3_kernel_v_read_readvariableop2
.savev2_adam_dense_3_bias_v_read_readvariableop4
0savev2_adam_dense_4_kernel_v_read_readvariableop2
.savev2_adam_dense_4_bias_v_read_readvariableop4
0savev2_adam_dense_5_kernel_v_read_readvariableop2
.savev2_adam_dense_5_bias_v_read_readvariableop
savev2_const

identity_1��MergeV2Checkpoints�
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Constl
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part2	
Const_1�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shard�
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename� 
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:2*
dtype0*�
value�B�2B0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUEB>layer_with_weights-0/optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB@layer_with_weights-0/optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB@layer_with_weights-0/optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-0/optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEBGlayer_with_weights-0/optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEBOlayer_with_weights-0/layer-0/layer-0/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEBOlayer_with_weights-0/layer-0/layer-1/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEBOlayer_with_weights-0/layer-0/layer-2/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEBOlayer_with_weights-0/layer-0/layer-3/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEBIlayer_with_weights-0/keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEBIlayer_with_weights-0/keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBIlayer_with_weights-0/keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEBIlayer_with_weights-0/keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/0/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/1/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/2/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/3/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/4/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/5/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/6/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/7/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/8/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/9/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBbtrainable_variables/10/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBbtrainable_variables/11/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/0/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/1/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/2/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/3/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/4/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/5/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/6/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/7/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/8/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/9/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBbtrainable_variables/10/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBbtrainable_variables/11/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_names�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:2*
dtype0*w
valuenBl2B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slices�
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0*savev2_conv2d_3_kernel_read_readvariableop(savev2_conv2d_3_bias_read_readvariableop*savev2_conv2d_4_kernel_read_readvariableop(savev2_conv2d_4_bias_read_readvariableop*savev2_conv2d_5_kernel_read_readvariableop(savev2_conv2d_5_bias_read_readvariableop)savev2_dense_3_kernel_read_readvariableop'savev2_dense_3_bias_read_readvariableop)savev2_dense_4_kernel_read_readvariableop'savev2_dense_4_bias_read_readvariableop)savev2_dense_5_kernel_read_readvariableop'savev2_dense_5_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop#savev2_variable_read_readvariableop%savev2_variable_1_read_readvariableop%savev2_variable_2_read_readvariableop%savev2_variable_3_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop1savev2_adam_conv2d_3_kernel_m_read_readvariableop/savev2_adam_conv2d_3_bias_m_read_readvariableop1savev2_adam_conv2d_4_kernel_m_read_readvariableop/savev2_adam_conv2d_4_bias_m_read_readvariableop1savev2_adam_conv2d_5_kernel_m_read_readvariableop/savev2_adam_conv2d_5_bias_m_read_readvariableop0savev2_adam_dense_3_kernel_m_read_readvariableop.savev2_adam_dense_3_bias_m_read_readvariableop0savev2_adam_dense_4_kernel_m_read_readvariableop.savev2_adam_dense_4_bias_m_read_readvariableop0savev2_adam_dense_5_kernel_m_read_readvariableop.savev2_adam_dense_5_bias_m_read_readvariableop1savev2_adam_conv2d_3_kernel_v_read_readvariableop/savev2_adam_conv2d_3_bias_v_read_readvariableop1savev2_adam_conv2d_4_kernel_v_read_readvariableop/savev2_adam_conv2d_4_bias_v_read_readvariableop1savev2_adam_conv2d_5_kernel_v_read_readvariableop/savev2_adam_conv2d_5_bias_v_read_readvariableop0savev2_adam_dense_3_kernel_v_read_readvariableop.savev2_adam_dense_3_bias_v_read_readvariableop0savev2_adam_dense_4_kernel_v_read_readvariableop.savev2_adam_dense_4_bias_v_read_readvariableop0savev2_adam_dense_5_kernel_v_read_readvariableop.savev2_adam_dense_5_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *@
dtypes6
422					2
SaveV2�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixes�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identitym

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*�
_input_shapes�
�: ::::: : :
��:�:	�:::: : : : : ::::: : : : ::::: : :
��:�:	�:::::::: : :
��:�:	�:::: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:,(
&
_output_shapes
:: 

_output_shapes
::,(
&
_output_shapes
:: 

_output_shapes
::,(
&
_output_shapes
: : 

_output_shapes
: :&"
 
_output_shapes
:
��:!

_output_shapes	
:�:%	!

_output_shapes
:	�: 


_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: : 

_output_shapes
:: 

_output_shapes
:: 

_output_shapes
:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :,(
&
_output_shapes
:: 

_output_shapes
::,(
&
_output_shapes
:: 

_output_shapes
::,(
&
_output_shapes
: : 

_output_shapes
: :& "
 
_output_shapes
:
��:!!

_output_shapes	
:�:%"!

_output_shapes
:	�: #

_output_shapes
::$$ 

_output_shapes

:: %

_output_shapes
::,&(
&
_output_shapes
:: '

_output_shapes
::,((
&
_output_shapes
:: )

_output_shapes
::,*(
&
_output_shapes
: : +

_output_shapes
: :&,"
 
_output_shapes
:
��:!-

_output_shapes	
:�:%.!

_output_shapes
:	�: /

_output_shapes
::$0 

_output_shapes

:: 1

_output_shapes
::2

_output_shapes
: 
�
b
D__inference_dropout_6_layer_call_and_return_conditional_losses_33985

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:����������2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:����������2

Identity_1"!

identity_1Identity_1:output:0*'
_input_shapes
:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
}
(__inference_conv2d_5_layer_call_fn_33925

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_5_layer_call_and_return_conditional_losses_315052
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*/
_output_shapes
:���������	 2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������	::22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������	
 
_user_specified_nameinputs
�

�
C__inference_conv2d_5_layer_call_and_return_conditional_losses_31505

inputs"
conv2d_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: *
dtype02
Conv2D/ReadVariableOp�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 *
paddingSAME*
strides
2
Conv2D�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:���������	 2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*
T0*/
_output_shapes
:���������	 2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������	::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:���������	
 
_user_specified_nameinputs
�S
�
G__inference_sequential_4_layer_call_and_return_conditional_losses_31761
sequential_3_input
conv2d_3_31710
conv2d_3_31712
conv2d_4_31717
conv2d_4_31719
conv2d_5_31724
conv2d_5_31726
dense_3_31731
dense_3_31733
dense_4_31737
dense_4_31739
dense_5_31743
dense_5_31745
identity�� conv2d_3/StatefulPartitionedCall� conv2d_4/StatefulPartitionedCall� conv2d_5/StatefulPartitionedCall�dense_3/StatefulPartitionedCall�0dense_3/kernel/Regularizer/Square/ReadVariableOp�dense_4/StatefulPartitionedCall�0dense_4/kernel/Regularizer/Square/ReadVariableOp�dense_5/StatefulPartitionedCall�
sequential_3/PartitionedCallPartitionedCallsequential_3_input*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_313082
sequential_3/PartitionedCallm
rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2
rescaling_3/Cast/xq
rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rescaling_3/Cast_1/x�
rescaling_3/mulMul%sequential_3/PartitionedCall:output:0rescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/mul�
rescaling_3/addAddV2rescaling_3/mul:z:0rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/add�
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCallrescaling_3/add:z:0conv2d_3_31710conv2d_3_31712*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_313892"
 conv2d_3/StatefulPartitionedCall�
dropout_4/PartitionedCallPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_4_layer_call_and_return_conditional_losses_314222
dropout_4/PartitionedCall�
max_pooling2d_3/PartitionedCallPartitionedCall"dropout_4/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_313172!
max_pooling2d_3/PartitionedCall�
 conv2d_4/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_3/PartitionedCall:output:0conv2d_4_31717conv2d_4_31719*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_314472"
 conv2d_4/StatefulPartitionedCall�
dropout_5/PartitionedCallPartitionedCall)conv2d_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_5_layer_call_and_return_conditional_losses_314802
dropout_5/PartitionedCall�
max_pooling2d_4/PartitionedCallPartitionedCall"dropout_5/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_4_layer_call_and_return_conditional_losses_313292!
max_pooling2d_4/PartitionedCall�
 conv2d_5/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_4/PartitionedCall:output:0conv2d_5_31724conv2d_5_31726*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_5_layer_call_and_return_conditional_losses_315052"
 conv2d_5/StatefulPartitionedCall�
max_pooling2d_5/PartitionedCallPartitionedCall)conv2d_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_5_layer_call_and_return_conditional_losses_313412!
max_pooling2d_5/PartitionedCall�
flatten_1/PartitionedCallPartitionedCall(max_pooling2d_5/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_flatten_1_layer_call_and_return_conditional_losses_315282
flatten_1/PartitionedCall�
dense_3/StatefulPartitionedCallStatefulPartitionedCall"flatten_1/PartitionedCall:output:0dense_3_31731dense_3_31733*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_315532!
dense_3/StatefulPartitionedCall�
dropout_6/PartitionedCallPartitionedCall(dense_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_6_layer_call_and_return_conditional_losses_315862
dropout_6/PartitionedCall�
dense_4/StatefulPartitionedCallStatefulPartitionedCall"dropout_6/PartitionedCall:output:0dense_4_31737dense_4_31739*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_4_layer_call_and_return_conditional_losses_316162!
dense_4/StatefulPartitionedCall�
dropout_7/PartitionedCallPartitionedCall(dense_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_7_layer_call_and_return_conditional_losses_316492
dropout_7/PartitionedCall�
dense_5/StatefulPartitionedCallStatefulPartitionedCall"dropout_7/PartitionedCall:output:0dense_5_31743dense_5_31745*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_5_layer_call_and_return_conditional_losses_316732!
dense_5/StatefulPartitionedCall�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_31731* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_4_31737*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity(dense_5/StatefulPartitionedCall:output:0!^conv2d_3/StatefulPartitionedCall!^conv2d_4/StatefulPartitionedCall!^conv2d_5/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall1^dense_3/kernel/Regularizer/Square/ReadVariableOp ^dense_4/StatefulPartitionedCall1^dense_4/kernel/Regularizer/Square/ReadVariableOp ^dense_5/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2D
 conv2d_4/StatefulPartitionedCall conv2d_4/StatefulPartitionedCall2D
 conv2d_5/StatefulPartitionedCall conv2d_5/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2B
dense_5/StatefulPartitionedCalldense_5/StatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_3_input
�
b
D__inference_dropout_4_layer_call_and_return_conditional_losses_31422

inputs

identity_1b
IdentityIdentityinputs*
T0*/
_output_shapes
:���������&2

Identityq

Identity_1IdentityIdentity:output:0*
T0*/
_output_shapes
:���������&2

Identity_1"!

identity_1Identity_1:output:0*.
_input_shapes
:���������&:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
`
D__inference_softmax_1_layer_call_and_return_conditional_losses_33458

inputs
identityW
SoftmaxSoftmaxinputs*
T0*'
_output_shapes
:���������2	
Softmaxe
IdentityIdentitySoftmax:softmax:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*&
_input_shapes
:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
}
(__inference_conv2d_3_layer_call_fn_33831

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_313892
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������&::22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�

�
,__inference_sequential_4_layer_call_fn_33424

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_4_layer_call_and_return_conditional_losses_318292
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
U
,__inference_sequential_3_layer_call_fn_31311
random_flip_1_input
identity�
PartitionedCallPartitionedCallrandom_flip_1_input*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_313082
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&:d `
/
_output_shapes
:���������&
-
_user_specified_namerandom_flip_1_input
�
c
D__inference_dropout_5_layer_call_and_return_conditional_losses_31475

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Const{
dropout/MulMulinputsdropout/Const:output:0*
T0*/
_output_shapes
:���������2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*/
_output_shapes
:���������*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*/
_output_shapes
:���������2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*/
_output_shapes
:���������2
dropout/Cast�
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*/
_output_shapes
:���������2
dropout/Mul_1m
IdentityIdentitydropout/Mul_1:z:0*
T0*/
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
�
}
(__inference_conv2d_4_layer_call_fn_33878

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_314472
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*/
_output_shapes
:���������2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������::22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
C__inference_conv2d_3_layer_call_and_return_conditional_losses_31389

inputs"
conv2d_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&*
paddingSAME*
strides
2
Conv2D�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:���������&2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������&::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
b
)__inference_dropout_5_layer_call_fn_33900

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_5_layer_call_and_return_conditional_losses_314752
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*/
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
�
b
)__inference_dropout_6_layer_call_fn_33990

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_6_layer_call_and_return_conditional_losses_315812
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*'
_input_shapes
:����������22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
b
D__inference_dropout_5_layer_call_and_return_conditional_losses_33895

inputs

identity_1b
IdentityIdentityinputs*
T0*/
_output_shapes
:���������2

Identityq

Identity_1IdentityIdentity:output:0*
T0*/
_output_shapes
:���������2

Identity_1"!

identity_1Identity_1:output:0*.
_input_shapes
:���������:W S
/
_output_shapes
:���������
 
_user_specified_nameinputs
�
|
'__inference_dense_5_layer_call_fn_34074

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_5_layer_call_and_return_conditional_losses_316732
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
b
D__inference_dropout_7_layer_call_and_return_conditional_losses_31649

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:���������2

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:���������2

Identity_1"!

identity_1Identity_1:output:0*&
_input_shapes
:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
̧
�
G__inference_sequential_4_layer_call_and_return_conditional_losses_33318

inputsL
Hsequential_3_random_rotation_1_stateful_uniform_statefuluniform_resourceO
Ksequential_3_random_translation_1_stateful_uniform_statefuluniform_resourceH
Dsequential_3_random_zoom_1_stateful_uniform_statefuluniform_resource+
'conv2d_3_conv2d_readvariableop_resource,
(conv2d_3_biasadd_readvariableop_resource+
'conv2d_4_conv2d_readvariableop_resource,
(conv2d_4_biasadd_readvariableop_resource+
'conv2d_5_conv2d_readvariableop_resource,
(conv2d_5_biasadd_readvariableop_resource*
&dense_3_matmul_readvariableop_resource+
'dense_3_biasadd_readvariableop_resource*
&dense_4_matmul_readvariableop_resource+
'dense_4_biasadd_readvariableop_resource*
&dense_5_matmul_readvariableop_resource+
'dense_5_biasadd_readvariableop_resource
identity��conv2d_3/BiasAdd/ReadVariableOp�conv2d_3/Conv2D/ReadVariableOp�conv2d_4/BiasAdd/ReadVariableOp�conv2d_4/Conv2D/ReadVariableOp�conv2d_5/BiasAdd/ReadVariableOp�conv2d_5/Conv2D/ReadVariableOp�dense_3/BiasAdd/ReadVariableOp�dense_3/MatMul/ReadVariableOp�0dense_3/kernel/Regularizer/Square/ReadVariableOp�dense_4/BiasAdd/ReadVariableOp�dense_4/MatMul/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�dense_5/BiasAdd/ReadVariableOp�dense_5/MatMul/ReadVariableOp�?sequential_3/random_rotation_1/stateful_uniform/StatefulUniform�Bsequential_3/random_translation_1/stateful_uniform/StatefulUniform�Dsequential_3/random_translation_1/stateful_uniform_1/StatefulUniform�;sequential_3/random_zoom_1/stateful_uniform/StatefulUniform�
Dsequential_3/random_flip_1/random_flip_left_right/control_dependencyIdentityinputs*
T0*
_class
loc:@inputs*/
_output_shapes
:���������&2F
Dsequential_3/random_flip_1/random_flip_left_right/control_dependency�
7sequential_3/random_flip_1/random_flip_left_right/ShapeShapeMsequential_3/random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:29
7sequential_3/random_flip_1/random_flip_left_right/Shape�
Esequential_3/random_flip_1/random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2G
Esequential_3/random_flip_1/random_flip_left_right/strided_slice/stack�
Gsequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2I
Gsequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_1�
Gsequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2I
Gsequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_2�
?sequential_3/random_flip_1/random_flip_left_right/strided_sliceStridedSlice@sequential_3/random_flip_1/random_flip_left_right/Shape:output:0Nsequential_3/random_flip_1/random_flip_left_right/strided_slice/stack:output:0Psequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_1:output:0Psequential_3/random_flip_1/random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2A
?sequential_3/random_flip_1/random_flip_left_right/strided_slice�
Fsequential_3/random_flip_1/random_flip_left_right/random_uniform/shapePackHsequential_3/random_flip_1/random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2H
Fsequential_3/random_flip_1/random_flip_left_right/random_uniform/shape�
Dsequential_3/random_flip_1/random_flip_left_right/random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2F
Dsequential_3/random_flip_1/random_flip_left_right/random_uniform/min�
Dsequential_3/random_flip_1/random_flip_left_right/random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2F
Dsequential_3/random_flip_1/random_flip_left_right/random_uniform/max�
Nsequential_3/random_flip_1/random_flip_left_right/random_uniform/RandomUniformRandomUniformOsequential_3/random_flip_1/random_flip_left_right/random_uniform/shape:output:0*
T0*#
_output_shapes
:���������*
dtype02P
Nsequential_3/random_flip_1/random_flip_left_right/random_uniform/RandomUniform�
Dsequential_3/random_flip_1/random_flip_left_right/random_uniform/MulMulWsequential_3/random_flip_1/random_flip_left_right/random_uniform/RandomUniform:output:0Msequential_3/random_flip_1/random_flip_left_right/random_uniform/max:output:0*
T0*#
_output_shapes
:���������2F
Dsequential_3/random_flip_1/random_flip_left_right/random_uniform/Mul�
Asequential_3/random_flip_1/random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2C
Asequential_3/random_flip_1/random_flip_left_right/Reshape/shape/1�
Asequential_3/random_flip_1/random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2C
Asequential_3/random_flip_1/random_flip_left_right/Reshape/shape/2�
Asequential_3/random_flip_1/random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2C
Asequential_3/random_flip_1/random_flip_left_right/Reshape/shape/3�
?sequential_3/random_flip_1/random_flip_left_right/Reshape/shapePackHsequential_3/random_flip_1/random_flip_left_right/strided_slice:output:0Jsequential_3/random_flip_1/random_flip_left_right/Reshape/shape/1:output:0Jsequential_3/random_flip_1/random_flip_left_right/Reshape/shape/2:output:0Jsequential_3/random_flip_1/random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2A
?sequential_3/random_flip_1/random_flip_left_right/Reshape/shape�
9sequential_3/random_flip_1/random_flip_left_right/ReshapeReshapeHsequential_3/random_flip_1/random_flip_left_right/random_uniform/Mul:z:0Hsequential_3/random_flip_1/random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:���������2;
9sequential_3/random_flip_1/random_flip_left_right/Reshape�
7sequential_3/random_flip_1/random_flip_left_right/RoundRoundBsequential_3/random_flip_1/random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:���������29
7sequential_3/random_flip_1/random_flip_left_right/Round�
@sequential_3/random_flip_1/random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:2B
@sequential_3/random_flip_1/random_flip_left_right/ReverseV2/axis�
;sequential_3/random_flip_1/random_flip_left_right/ReverseV2	ReverseV2Msequential_3/random_flip_1/random_flip_left_right/control_dependency:output:0Isequential_3/random_flip_1/random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:���������&2=
;sequential_3/random_flip_1/random_flip_left_right/ReverseV2�
5sequential_3/random_flip_1/random_flip_left_right/mulMul;sequential_3/random_flip_1/random_flip_left_right/Round:y:0Dsequential_3/random_flip_1/random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:���������&27
5sequential_3/random_flip_1/random_flip_left_right/mul�
7sequential_3/random_flip_1/random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?29
7sequential_3/random_flip_1/random_flip_left_right/sub/x�
5sequential_3/random_flip_1/random_flip_left_right/subSub@sequential_3/random_flip_1/random_flip_left_right/sub/x:output:0;sequential_3/random_flip_1/random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:���������27
5sequential_3/random_flip_1/random_flip_left_right/sub�
7sequential_3/random_flip_1/random_flip_left_right/mul_1Mul9sequential_3/random_flip_1/random_flip_left_right/sub:z:0Msequential_3/random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:���������&29
7sequential_3/random_flip_1/random_flip_left_right/mul_1�
5sequential_3/random_flip_1/random_flip_left_right/addAddV29sequential_3/random_flip_1/random_flip_left_right/mul:z:0;sequential_3/random_flip_1/random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:���������&27
5sequential_3/random_flip_1/random_flip_left_right/add�
$sequential_3/random_rotation_1/ShapeShape9sequential_3/random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:2&
$sequential_3/random_rotation_1/Shape�
2sequential_3/random_rotation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 24
2sequential_3/random_rotation_1/strided_slice/stack�
4sequential_3/random_rotation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:26
4sequential_3/random_rotation_1/strided_slice/stack_1�
4sequential_3/random_rotation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:26
4sequential_3/random_rotation_1/strided_slice/stack_2�
,sequential_3/random_rotation_1/strided_sliceStridedSlice-sequential_3/random_rotation_1/Shape:output:0;sequential_3/random_rotation_1/strided_slice/stack:output:0=sequential_3/random_rotation_1/strided_slice/stack_1:output:0=sequential_3/random_rotation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2.
,sequential_3/random_rotation_1/strided_slice�
4sequential_3/random_rotation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:26
4sequential_3/random_rotation_1/strided_slice_1/stack�
6sequential_3/random_rotation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:28
6sequential_3/random_rotation_1/strided_slice_1/stack_1�
6sequential_3/random_rotation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:28
6sequential_3/random_rotation_1/strided_slice_1/stack_2�
.sequential_3/random_rotation_1/strided_slice_1StridedSlice-sequential_3/random_rotation_1/Shape:output:0=sequential_3/random_rotation_1/strided_slice_1/stack:output:0?sequential_3/random_rotation_1/strided_slice_1/stack_1:output:0?sequential_3/random_rotation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask20
.sequential_3/random_rotation_1/strided_slice_1�
#sequential_3/random_rotation_1/CastCast7sequential_3/random_rotation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2%
#sequential_3/random_rotation_1/Cast�
4sequential_3/random_rotation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:26
4sequential_3/random_rotation_1/strided_slice_2/stack�
6sequential_3/random_rotation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:28
6sequential_3/random_rotation_1/strided_slice_2/stack_1�
6sequential_3/random_rotation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:28
6sequential_3/random_rotation_1/strided_slice_2/stack_2�
.sequential_3/random_rotation_1/strided_slice_2StridedSlice-sequential_3/random_rotation_1/Shape:output:0=sequential_3/random_rotation_1/strided_slice_2/stack:output:0?sequential_3/random_rotation_1/strided_slice_2/stack_1:output:0?sequential_3/random_rotation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask20
.sequential_3/random_rotation_1/strided_slice_2�
%sequential_3/random_rotation_1/Cast_1Cast7sequential_3/random_rotation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2'
%sequential_3/random_rotation_1/Cast_1�
5sequential_3/random_rotation_1/stateful_uniform/shapePack5sequential_3/random_rotation_1/strided_slice:output:0*
N*
T0*
_output_shapes
:27
5sequential_3/random_rotation_1/stateful_uniform/shape�
3sequential_3/random_rotation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *   �25
3sequential_3/random_rotation_1/stateful_uniform/min�
3sequential_3/random_rotation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    25
3sequential_3/random_rotation_1/stateful_uniform/max�
Isequential_3/random_rotation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2K
Isequential_3/random_rotation_1/stateful_uniform/StatefulUniform/algorithm�
?sequential_3/random_rotation_1/stateful_uniform/StatefulUniformStatefulUniformHsequential_3_random_rotation_1_stateful_uniform_statefuluniform_resourceRsequential_3/random_rotation_1/stateful_uniform/StatefulUniform/algorithm:output:0>sequential_3/random_rotation_1/stateful_uniform/shape:output:0*#
_output_shapes
:���������*
shape_dtype02A
?sequential_3/random_rotation_1/stateful_uniform/StatefulUniform�
3sequential_3/random_rotation_1/stateful_uniform/subSub<sequential_3/random_rotation_1/stateful_uniform/max:output:0<sequential_3/random_rotation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 25
3sequential_3/random_rotation_1/stateful_uniform/sub�
3sequential_3/random_rotation_1/stateful_uniform/mulMulHsequential_3/random_rotation_1/stateful_uniform/StatefulUniform:output:07sequential_3/random_rotation_1/stateful_uniform/sub:z:0*
T0*#
_output_shapes
:���������25
3sequential_3/random_rotation_1/stateful_uniform/mul�
/sequential_3/random_rotation_1/stateful_uniformAdd7sequential_3/random_rotation_1/stateful_uniform/mul:z:0<sequential_3/random_rotation_1/stateful_uniform/min:output:0*
T0*#
_output_shapes
:���������21
/sequential_3/random_rotation_1/stateful_uniform�
4sequential_3/random_rotation_1/rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?26
4sequential_3/random_rotation_1/rotation_matrix/sub/y�
2sequential_3/random_rotation_1/rotation_matrix/subSub)sequential_3/random_rotation_1/Cast_1:y:0=sequential_3/random_rotation_1/rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 24
2sequential_3/random_rotation_1/rotation_matrix/sub�
2sequential_3/random_rotation_1/rotation_matrix/CosCos3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������24
2sequential_3/random_rotation_1/rotation_matrix/Cos�
6sequential_3/random_rotation_1/rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?28
6sequential_3/random_rotation_1/rotation_matrix/sub_1/y�
4sequential_3/random_rotation_1/rotation_matrix/sub_1Sub)sequential_3/random_rotation_1/Cast_1:y:0?sequential_3/random_rotation_1/rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 26
4sequential_3/random_rotation_1/rotation_matrix/sub_1�
2sequential_3/random_rotation_1/rotation_matrix/mulMul6sequential_3/random_rotation_1/rotation_matrix/Cos:y:08sequential_3/random_rotation_1/rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:���������24
2sequential_3/random_rotation_1/rotation_matrix/mul�
2sequential_3/random_rotation_1/rotation_matrix/SinSin3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������24
2sequential_3/random_rotation_1/rotation_matrix/Sin�
6sequential_3/random_rotation_1/rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?28
6sequential_3/random_rotation_1/rotation_matrix/sub_2/y�
4sequential_3/random_rotation_1/rotation_matrix/sub_2Sub'sequential_3/random_rotation_1/Cast:y:0?sequential_3/random_rotation_1/rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 26
4sequential_3/random_rotation_1/rotation_matrix/sub_2�
4sequential_3/random_rotation_1/rotation_matrix/mul_1Mul6sequential_3/random_rotation_1/rotation_matrix/Sin:y:08sequential_3/random_rotation_1/rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/mul_1�
4sequential_3/random_rotation_1/rotation_matrix/sub_3Sub6sequential_3/random_rotation_1/rotation_matrix/mul:z:08sequential_3/random_rotation_1/rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/sub_3�
4sequential_3/random_rotation_1/rotation_matrix/sub_4Sub6sequential_3/random_rotation_1/rotation_matrix/sub:z:08sequential_3/random_rotation_1/rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/sub_4�
8sequential_3/random_rotation_1/rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2:
8sequential_3/random_rotation_1/rotation_matrix/truediv/y�
6sequential_3/random_rotation_1/rotation_matrix/truedivRealDiv8sequential_3/random_rotation_1/rotation_matrix/sub_4:z:0Asequential_3/random_rotation_1/rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:���������28
6sequential_3/random_rotation_1/rotation_matrix/truediv�
6sequential_3/random_rotation_1/rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?28
6sequential_3/random_rotation_1/rotation_matrix/sub_5/y�
4sequential_3/random_rotation_1/rotation_matrix/sub_5Sub'sequential_3/random_rotation_1/Cast:y:0?sequential_3/random_rotation_1/rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 26
4sequential_3/random_rotation_1/rotation_matrix/sub_5�
4sequential_3/random_rotation_1/rotation_matrix/Sin_1Sin3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/Sin_1�
6sequential_3/random_rotation_1/rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?28
6sequential_3/random_rotation_1/rotation_matrix/sub_6/y�
4sequential_3/random_rotation_1/rotation_matrix/sub_6Sub)sequential_3/random_rotation_1/Cast_1:y:0?sequential_3/random_rotation_1/rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 26
4sequential_3/random_rotation_1/rotation_matrix/sub_6�
4sequential_3/random_rotation_1/rotation_matrix/mul_2Mul8sequential_3/random_rotation_1/rotation_matrix/Sin_1:y:08sequential_3/random_rotation_1/rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/mul_2�
4sequential_3/random_rotation_1/rotation_matrix/Cos_1Cos3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/Cos_1�
6sequential_3/random_rotation_1/rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?28
6sequential_3/random_rotation_1/rotation_matrix/sub_7/y�
4sequential_3/random_rotation_1/rotation_matrix/sub_7Sub'sequential_3/random_rotation_1/Cast:y:0?sequential_3/random_rotation_1/rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 26
4sequential_3/random_rotation_1/rotation_matrix/sub_7�
4sequential_3/random_rotation_1/rotation_matrix/mul_3Mul8sequential_3/random_rotation_1/rotation_matrix/Cos_1:y:08sequential_3/random_rotation_1/rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/mul_3�
2sequential_3/random_rotation_1/rotation_matrix/addAddV28sequential_3/random_rotation_1/rotation_matrix/mul_2:z:08sequential_3/random_rotation_1/rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:���������24
2sequential_3/random_rotation_1/rotation_matrix/add�
4sequential_3/random_rotation_1/rotation_matrix/sub_8Sub8sequential_3/random_rotation_1/rotation_matrix/sub_5:z:06sequential_3/random_rotation_1/rotation_matrix/add:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/sub_8�
:sequential_3/random_rotation_1/rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2<
:sequential_3/random_rotation_1/rotation_matrix/truediv_1/y�
8sequential_3/random_rotation_1/rotation_matrix/truediv_1RealDiv8sequential_3/random_rotation_1/rotation_matrix/sub_8:z:0Csequential_3/random_rotation_1/rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:���������2:
8sequential_3/random_rotation_1/rotation_matrix/truediv_1�
4sequential_3/random_rotation_1/rotation_matrix/ShapeShape3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*
_output_shapes
:26
4sequential_3/random_rotation_1/rotation_matrix/Shape�
Bsequential_3/random_rotation_1/rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2D
Bsequential_3/random_rotation_1/rotation_matrix/strided_slice/stack�
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2F
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_1�
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2F
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_2�
<sequential_3/random_rotation_1/rotation_matrix/strided_sliceStridedSlice=sequential_3/random_rotation_1/rotation_matrix/Shape:output:0Ksequential_3/random_rotation_1/rotation_matrix/strided_slice/stack:output:0Msequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_1:output:0Msequential_3/random_rotation_1/rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2>
<sequential_3/random_rotation_1/rotation_matrix/strided_slice�
4sequential_3/random_rotation_1/rotation_matrix/Cos_2Cos3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/Cos_2�
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        2F
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_1�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_2�
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_1StridedSlice8sequential_3/random_rotation_1/rotation_matrix/Cos_2:y:0Msequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_1:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2@
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_1�
4sequential_3/random_rotation_1/rotation_matrix/Sin_2Sin3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/Sin_2�
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        2F
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_1�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_2�
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_2StridedSlice8sequential_3/random_rotation_1/rotation_matrix/Sin_2:y:0Msequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_1:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2@
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_2�
2sequential_3/random_rotation_1/rotation_matrix/NegNegGsequential_3/random_rotation_1/rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������24
2sequential_3/random_rotation_1/rotation_matrix/Neg�
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        2F
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_1�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_2�
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_3StridedSlice:sequential_3/random_rotation_1/rotation_matrix/truediv:z:0Msequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_1:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2@
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_3�
4sequential_3/random_rotation_1/rotation_matrix/Sin_3Sin3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/Sin_3�
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        2F
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_1�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_2�
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_4StridedSlice8sequential_3/random_rotation_1/rotation_matrix/Sin_3:y:0Msequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_1:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2@
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_4�
4sequential_3/random_rotation_1/rotation_matrix/Cos_3Cos3sequential_3/random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/Cos_3�
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        2F
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_1�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_2�
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_5StridedSlice8sequential_3/random_rotation_1/rotation_matrix/Cos_3:y:0Msequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_1:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2@
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_5�
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        2F
Dsequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_1�
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2H
Fsequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_2�
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_6StridedSlice<sequential_3/random_rotation_1/rotation_matrix/truediv_1:z:0Msequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_1:output:0Osequential_3/random_rotation_1/rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask2@
>sequential_3/random_rotation_1/rotation_matrix/strided_slice_6�
:sequential_3/random_rotation_1/rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2<
:sequential_3/random_rotation_1/rotation_matrix/zeros/mul/y�
8sequential_3/random_rotation_1/rotation_matrix/zeros/mulMulEsequential_3/random_rotation_1/rotation_matrix/strided_slice:output:0Csequential_3/random_rotation_1/rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2:
8sequential_3/random_rotation_1/rotation_matrix/zeros/mul�
;sequential_3/random_rotation_1/rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2=
;sequential_3/random_rotation_1/rotation_matrix/zeros/Less/y�
9sequential_3/random_rotation_1/rotation_matrix/zeros/LessLess<sequential_3/random_rotation_1/rotation_matrix/zeros/mul:z:0Dsequential_3/random_rotation_1/rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2;
9sequential_3/random_rotation_1/rotation_matrix/zeros/Less�
=sequential_3/random_rotation_1/rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2?
=sequential_3/random_rotation_1/rotation_matrix/zeros/packed/1�
;sequential_3/random_rotation_1/rotation_matrix/zeros/packedPackEsequential_3/random_rotation_1/rotation_matrix/strided_slice:output:0Fsequential_3/random_rotation_1/rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2=
;sequential_3/random_rotation_1/rotation_matrix/zeros/packed�
:sequential_3/random_rotation_1/rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2<
:sequential_3/random_rotation_1/rotation_matrix/zeros/Const�
4sequential_3/random_rotation_1/rotation_matrix/zerosFillDsequential_3/random_rotation_1/rotation_matrix/zeros/packed:output:0Csequential_3/random_rotation_1/rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������26
4sequential_3/random_rotation_1/rotation_matrix/zeros�
:sequential_3/random_rotation_1/rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2<
:sequential_3/random_rotation_1/rotation_matrix/concat/axis�
5sequential_3/random_rotation_1/rotation_matrix/concatConcatV2Gsequential_3/random_rotation_1/rotation_matrix/strided_slice_1:output:06sequential_3/random_rotation_1/rotation_matrix/Neg:y:0Gsequential_3/random_rotation_1/rotation_matrix/strided_slice_3:output:0Gsequential_3/random_rotation_1/rotation_matrix/strided_slice_4:output:0Gsequential_3/random_rotation_1/rotation_matrix/strided_slice_5:output:0Gsequential_3/random_rotation_1/rotation_matrix/strided_slice_6:output:0=sequential_3/random_rotation_1/rotation_matrix/zeros:output:0Csequential_3/random_rotation_1/rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������27
5sequential_3/random_rotation_1/rotation_matrix/concat�
.sequential_3/random_rotation_1/transform/ShapeShape9sequential_3/random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:20
.sequential_3/random_rotation_1/transform/Shape�
<sequential_3/random_rotation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2>
<sequential_3/random_rotation_1/transform/strided_slice/stack�
>sequential_3/random_rotation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2@
>sequential_3/random_rotation_1/transform/strided_slice/stack_1�
>sequential_3/random_rotation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2@
>sequential_3/random_rotation_1/transform/strided_slice/stack_2�
6sequential_3/random_rotation_1/transform/strided_sliceStridedSlice7sequential_3/random_rotation_1/transform/Shape:output:0Esequential_3/random_rotation_1/transform/strided_slice/stack:output:0Gsequential_3/random_rotation_1/transform/strided_slice/stack_1:output:0Gsequential_3/random_rotation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:28
6sequential_3/random_rotation_1/transform/strided_slice�
3sequential_3/random_rotation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    25
3sequential_3/random_rotation_1/transform/fill_value�
Csequential_3/random_rotation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV39sequential_3/random_flip_1/random_flip_left_right/add:z:0>sequential_3/random_rotation_1/rotation_matrix/concat:output:0?sequential_3/random_rotation_1/transform/strided_slice:output:0<sequential_3/random_rotation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2E
Csequential_3/random_rotation_1/transform/ImageProjectiveTransformV3�
'sequential_3/random_translation_1/ShapeShapeXsequential_3/random_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2)
'sequential_3/random_translation_1/Shape�
5sequential_3/random_translation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 27
5sequential_3/random_translation_1/strided_slice/stack�
7sequential_3/random_translation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:29
7sequential_3/random_translation_1/strided_slice/stack_1�
7sequential_3/random_translation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:29
7sequential_3/random_translation_1/strided_slice/stack_2�
/sequential_3/random_translation_1/strided_sliceStridedSlice0sequential_3/random_translation_1/Shape:output:0>sequential_3/random_translation_1/strided_slice/stack:output:0@sequential_3/random_translation_1/strided_slice/stack_1:output:0@sequential_3/random_translation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask21
/sequential_3/random_translation_1/strided_slice�
7sequential_3/random_translation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:29
7sequential_3/random_translation_1/strided_slice_1/stack�
9sequential_3/random_translation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2;
9sequential_3/random_translation_1/strided_slice_1/stack_1�
9sequential_3/random_translation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2;
9sequential_3/random_translation_1/strided_slice_1/stack_2�
1sequential_3/random_translation_1/strided_slice_1StridedSlice0sequential_3/random_translation_1/Shape:output:0@sequential_3/random_translation_1/strided_slice_1/stack:output:0Bsequential_3/random_translation_1/strided_slice_1/stack_1:output:0Bsequential_3/random_translation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask23
1sequential_3/random_translation_1/strided_slice_1�
&sequential_3/random_translation_1/CastCast:sequential_3/random_translation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2(
&sequential_3/random_translation_1/Cast�
7sequential_3/random_translation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:29
7sequential_3/random_translation_1/strided_slice_2/stack�
9sequential_3/random_translation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2;
9sequential_3/random_translation_1/strided_slice_2/stack_1�
9sequential_3/random_translation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2;
9sequential_3/random_translation_1/strided_slice_2/stack_2�
1sequential_3/random_translation_1/strided_slice_2StridedSlice0sequential_3/random_translation_1/Shape:output:0@sequential_3/random_translation_1/strided_slice_2/stack:output:0Bsequential_3/random_translation_1/strided_slice_2/stack_1:output:0Bsequential_3/random_translation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask23
1sequential_3/random_translation_1/strided_slice_2�
(sequential_3/random_translation_1/Cast_1Cast:sequential_3/random_translation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2*
(sequential_3/random_translation_1/Cast_1�
:sequential_3/random_translation_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2<
:sequential_3/random_translation_1/stateful_uniform/shape/1�
8sequential_3/random_translation_1/stateful_uniform/shapePack8sequential_3/random_translation_1/strided_slice:output:0Csequential_3/random_translation_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2:
8sequential_3/random_translation_1/stateful_uniform/shape�
6sequential_3/random_translation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    28
6sequential_3/random_translation_1/stateful_uniform/min�
6sequential_3/random_translation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *��>28
6sequential_3/random_translation_1/stateful_uniform/max�
Lsequential_3/random_translation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2N
Lsequential_3/random_translation_1/stateful_uniform/StatefulUniform/algorithm�
Bsequential_3/random_translation_1/stateful_uniform/StatefulUniformStatefulUniformKsequential_3_random_translation_1_stateful_uniform_statefuluniform_resourceUsequential_3/random_translation_1/stateful_uniform/StatefulUniform/algorithm:output:0Asequential_3/random_translation_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype02D
Bsequential_3/random_translation_1/stateful_uniform/StatefulUniform�
6sequential_3/random_translation_1/stateful_uniform/subSub?sequential_3/random_translation_1/stateful_uniform/max:output:0?sequential_3/random_translation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 28
6sequential_3/random_translation_1/stateful_uniform/sub�
6sequential_3/random_translation_1/stateful_uniform/mulMulKsequential_3/random_translation_1/stateful_uniform/StatefulUniform:output:0:sequential_3/random_translation_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������28
6sequential_3/random_translation_1/stateful_uniform/mul�
2sequential_3/random_translation_1/stateful_uniformAdd:sequential_3/random_translation_1/stateful_uniform/mul:z:0?sequential_3/random_translation_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������24
2sequential_3/random_translation_1/stateful_uniform�
%sequential_3/random_translation_1/mulMul6sequential_3/random_translation_1/stateful_uniform:z:0*sequential_3/random_translation_1/Cast:y:0*
T0*'
_output_shapes
:���������2'
%sequential_3/random_translation_1/mul�
<sequential_3/random_translation_1/stateful_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2>
<sequential_3/random_translation_1/stateful_uniform_1/shape/1�
:sequential_3/random_translation_1/stateful_uniform_1/shapePack8sequential_3/random_translation_1/strided_slice:output:0Esequential_3/random_translation_1/stateful_uniform_1/shape/1:output:0*
N*
T0*
_output_shapes
:2<
:sequential_3/random_translation_1/stateful_uniform_1/shape�
8sequential_3/random_translation_1/stateful_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2:
8sequential_3/random_translation_1/stateful_uniform_1/min�
8sequential_3/random_translation_1/stateful_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2:
8sequential_3/random_translation_1/stateful_uniform_1/max�
Nsequential_3/random_translation_1/stateful_uniform_1/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2P
Nsequential_3/random_translation_1/stateful_uniform_1/StatefulUniform/algorithm�
Dsequential_3/random_translation_1/stateful_uniform_1/StatefulUniformStatefulUniformKsequential_3_random_translation_1_stateful_uniform_statefuluniform_resourceWsequential_3/random_translation_1/stateful_uniform_1/StatefulUniform/algorithm:output:0Csequential_3/random_translation_1/stateful_uniform_1/shape:output:0C^sequential_3/random_translation_1/stateful_uniform/StatefulUniform*'
_output_shapes
:���������*
shape_dtype02F
Dsequential_3/random_translation_1/stateful_uniform_1/StatefulUniform�
8sequential_3/random_translation_1/stateful_uniform_1/subSubAsequential_3/random_translation_1/stateful_uniform_1/max:output:0Asequential_3/random_translation_1/stateful_uniform_1/min:output:0*
T0*
_output_shapes
: 2:
8sequential_3/random_translation_1/stateful_uniform_1/sub�
8sequential_3/random_translation_1/stateful_uniform_1/mulMulMsequential_3/random_translation_1/stateful_uniform_1/StatefulUniform:output:0<sequential_3/random_translation_1/stateful_uniform_1/sub:z:0*
T0*'
_output_shapes
:���������2:
8sequential_3/random_translation_1/stateful_uniform_1/mul�
4sequential_3/random_translation_1/stateful_uniform_1Add<sequential_3/random_translation_1/stateful_uniform_1/mul:z:0Asequential_3/random_translation_1/stateful_uniform_1/min:output:0*
T0*'
_output_shapes
:���������26
4sequential_3/random_translation_1/stateful_uniform_1�
'sequential_3/random_translation_1/mul_1Mul8sequential_3/random_translation_1/stateful_uniform_1:z:0,sequential_3/random_translation_1/Cast_1:y:0*
T0*'
_output_shapes
:���������2)
'sequential_3/random_translation_1/mul_1�
-sequential_3/random_translation_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2/
-sequential_3/random_translation_1/concat/axis�
(sequential_3/random_translation_1/concatConcatV2+sequential_3/random_translation_1/mul_1:z:0)sequential_3/random_translation_1/mul:z:06sequential_3/random_translation_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2*
(sequential_3/random_translation_1/concat�
:sequential_3/random_translation_1/translation_matrix/ShapeShape1sequential_3/random_translation_1/concat:output:0*
T0*
_output_shapes
:2<
:sequential_3/random_translation_1/translation_matrix/Shape�
Hsequential_3/random_translation_1/translation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2J
Hsequential_3/random_translation_1/translation_matrix/strided_slice/stack�
Jsequential_3/random_translation_1/translation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2L
Jsequential_3/random_translation_1/translation_matrix/strided_slice/stack_1�
Jsequential_3/random_translation_1/translation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2L
Jsequential_3/random_translation_1/translation_matrix/strided_slice/stack_2�
Bsequential_3/random_translation_1/translation_matrix/strided_sliceStridedSliceCsequential_3/random_translation_1/translation_matrix/Shape:output:0Qsequential_3/random_translation_1/translation_matrix/strided_slice/stack:output:0Ssequential_3/random_translation_1/translation_matrix/strided_slice/stack_1:output:0Ssequential_3/random_translation_1/translation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2D
Bsequential_3/random_translation_1/translation_matrix/strided_slice�
?sequential_3/random_translation_1/translation_matrix/ones/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2A
?sequential_3/random_translation_1/translation_matrix/ones/mul/y�
=sequential_3/random_translation_1/translation_matrix/ones/mulMulKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Hsequential_3/random_translation_1/translation_matrix/ones/mul/y:output:0*
T0*
_output_shapes
: 2?
=sequential_3/random_translation_1/translation_matrix/ones/mul�
@sequential_3/random_translation_1/translation_matrix/ones/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2B
@sequential_3/random_translation_1/translation_matrix/ones/Less/y�
>sequential_3/random_translation_1/translation_matrix/ones/LessLessAsequential_3/random_translation_1/translation_matrix/ones/mul:z:0Isequential_3/random_translation_1/translation_matrix/ones/Less/y:output:0*
T0*
_output_shapes
: 2@
>sequential_3/random_translation_1/translation_matrix/ones/Less�
Bsequential_3/random_translation_1/translation_matrix/ones/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2D
Bsequential_3/random_translation_1/translation_matrix/ones/packed/1�
@sequential_3/random_translation_1/translation_matrix/ones/packedPackKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Ksequential_3/random_translation_1/translation_matrix/ones/packed/1:output:0*
N*
T0*
_output_shapes
:2B
@sequential_3/random_translation_1/translation_matrix/ones/packed�
?sequential_3/random_translation_1/translation_matrix/ones/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2A
?sequential_3/random_translation_1/translation_matrix/ones/Const�
9sequential_3/random_translation_1/translation_matrix/onesFillIsequential_3/random_translation_1/translation_matrix/ones/packed:output:0Hsequential_3/random_translation_1/translation_matrix/ones/Const:output:0*
T0*'
_output_shapes
:���������2;
9sequential_3/random_translation_1/translation_matrix/ones�
@sequential_3/random_translation_1/translation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2B
@sequential_3/random_translation_1/translation_matrix/zeros/mul/y�
>sequential_3/random_translation_1/translation_matrix/zeros/mulMulKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Isequential_3/random_translation_1/translation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2@
>sequential_3/random_translation_1/translation_matrix/zeros/mul�
Asequential_3/random_translation_1/translation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2C
Asequential_3/random_translation_1/translation_matrix/zeros/Less/y�
?sequential_3/random_translation_1/translation_matrix/zeros/LessLessBsequential_3/random_translation_1/translation_matrix/zeros/mul:z:0Jsequential_3/random_translation_1/translation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2A
?sequential_3/random_translation_1/translation_matrix/zeros/Less�
Csequential_3/random_translation_1/translation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2E
Csequential_3/random_translation_1/translation_matrix/zeros/packed/1�
Asequential_3/random_translation_1/translation_matrix/zeros/packedPackKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Lsequential_3/random_translation_1/translation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2C
Asequential_3/random_translation_1/translation_matrix/zeros/packed�
@sequential_3/random_translation_1/translation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2B
@sequential_3/random_translation_1/translation_matrix/zeros/Const�
:sequential_3/random_translation_1/translation_matrix/zerosFillJsequential_3/random_translation_1/translation_matrix/zeros/packed:output:0Isequential_3/random_translation_1/translation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2<
:sequential_3/random_translation_1/translation_matrix/zeros�
Jsequential_3/random_translation_1/translation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2L
Jsequential_3/random_translation_1/translation_matrix/strided_slice_1/stack�
Lsequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2N
Lsequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_1�
Lsequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2N
Lsequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_2�
Dsequential_3/random_translation_1/translation_matrix/strided_slice_1StridedSlice1sequential_3/random_translation_1/concat:output:0Ssequential_3/random_translation_1/translation_matrix/strided_slice_1/stack:output:0Usequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_1:output:0Usequential_3/random_translation_1/translation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2F
Dsequential_3/random_translation_1/translation_matrix/strided_slice_1�
8sequential_3/random_translation_1/translation_matrix/NegNegMsequential_3/random_translation_1/translation_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2:
8sequential_3/random_translation_1/translation_matrix/Neg�
Bsequential_3/random_translation_1/translation_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2D
Bsequential_3/random_translation_1/translation_matrix/zeros_1/mul/y�
@sequential_3/random_translation_1/translation_matrix/zeros_1/mulMulKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Ksequential_3/random_translation_1/translation_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2B
@sequential_3/random_translation_1/translation_matrix/zeros_1/mul�
Csequential_3/random_translation_1/translation_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2E
Csequential_3/random_translation_1/translation_matrix/zeros_1/Less/y�
Asequential_3/random_translation_1/translation_matrix/zeros_1/LessLessDsequential_3/random_translation_1/translation_matrix/zeros_1/mul:z:0Lsequential_3/random_translation_1/translation_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2C
Asequential_3/random_translation_1/translation_matrix/zeros_1/Less�
Esequential_3/random_translation_1/translation_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2G
Esequential_3/random_translation_1/translation_matrix/zeros_1/packed/1�
Csequential_3/random_translation_1/translation_matrix/zeros_1/packedPackKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Nsequential_3/random_translation_1/translation_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2E
Csequential_3/random_translation_1/translation_matrix/zeros_1/packed�
Bsequential_3/random_translation_1/translation_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2D
Bsequential_3/random_translation_1/translation_matrix/zeros_1/Const�
<sequential_3/random_translation_1/translation_matrix/zeros_1FillLsequential_3/random_translation_1/translation_matrix/zeros_1/packed:output:0Ksequential_3/random_translation_1/translation_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������2>
<sequential_3/random_translation_1/translation_matrix/zeros_1�
Asequential_3/random_translation_1/translation_matrix/ones_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2C
Asequential_3/random_translation_1/translation_matrix/ones_1/mul/y�
?sequential_3/random_translation_1/translation_matrix/ones_1/mulMulKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Jsequential_3/random_translation_1/translation_matrix/ones_1/mul/y:output:0*
T0*
_output_shapes
: 2A
?sequential_3/random_translation_1/translation_matrix/ones_1/mul�
Bsequential_3/random_translation_1/translation_matrix/ones_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2D
Bsequential_3/random_translation_1/translation_matrix/ones_1/Less/y�
@sequential_3/random_translation_1/translation_matrix/ones_1/LessLessCsequential_3/random_translation_1/translation_matrix/ones_1/mul:z:0Ksequential_3/random_translation_1/translation_matrix/ones_1/Less/y:output:0*
T0*
_output_shapes
: 2B
@sequential_3/random_translation_1/translation_matrix/ones_1/Less�
Dsequential_3/random_translation_1/translation_matrix/ones_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2F
Dsequential_3/random_translation_1/translation_matrix/ones_1/packed/1�
Bsequential_3/random_translation_1/translation_matrix/ones_1/packedPackKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Msequential_3/random_translation_1/translation_matrix/ones_1/packed/1:output:0*
N*
T0*
_output_shapes
:2D
Bsequential_3/random_translation_1/translation_matrix/ones_1/packed�
Asequential_3/random_translation_1/translation_matrix/ones_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2C
Asequential_3/random_translation_1/translation_matrix/ones_1/Const�
;sequential_3/random_translation_1/translation_matrix/ones_1FillKsequential_3/random_translation_1/translation_matrix/ones_1/packed:output:0Jsequential_3/random_translation_1/translation_matrix/ones_1/Const:output:0*
T0*'
_output_shapes
:���������2=
;sequential_3/random_translation_1/translation_matrix/ones_1�
Jsequential_3/random_translation_1/translation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2L
Jsequential_3/random_translation_1/translation_matrix/strided_slice_2/stack�
Lsequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2N
Lsequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_1�
Lsequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2N
Lsequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_2�
Dsequential_3/random_translation_1/translation_matrix/strided_slice_2StridedSlice1sequential_3/random_translation_1/concat:output:0Ssequential_3/random_translation_1/translation_matrix/strided_slice_2/stack:output:0Usequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_1:output:0Usequential_3/random_translation_1/translation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2F
Dsequential_3/random_translation_1/translation_matrix/strided_slice_2�
:sequential_3/random_translation_1/translation_matrix/Neg_1NegMsequential_3/random_translation_1/translation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2<
:sequential_3/random_translation_1/translation_matrix/Neg_1�
Bsequential_3/random_translation_1/translation_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2D
Bsequential_3/random_translation_1/translation_matrix/zeros_2/mul/y�
@sequential_3/random_translation_1/translation_matrix/zeros_2/mulMulKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Ksequential_3/random_translation_1/translation_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2B
@sequential_3/random_translation_1/translation_matrix/zeros_2/mul�
Csequential_3/random_translation_1/translation_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2E
Csequential_3/random_translation_1/translation_matrix/zeros_2/Less/y�
Asequential_3/random_translation_1/translation_matrix/zeros_2/LessLessDsequential_3/random_translation_1/translation_matrix/zeros_2/mul:z:0Lsequential_3/random_translation_1/translation_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2C
Asequential_3/random_translation_1/translation_matrix/zeros_2/Less�
Esequential_3/random_translation_1/translation_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2G
Esequential_3/random_translation_1/translation_matrix/zeros_2/packed/1�
Csequential_3/random_translation_1/translation_matrix/zeros_2/packedPackKsequential_3/random_translation_1/translation_matrix/strided_slice:output:0Nsequential_3/random_translation_1/translation_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2E
Csequential_3/random_translation_1/translation_matrix/zeros_2/packed�
Bsequential_3/random_translation_1/translation_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2D
Bsequential_3/random_translation_1/translation_matrix/zeros_2/Const�
<sequential_3/random_translation_1/translation_matrix/zeros_2FillLsequential_3/random_translation_1/translation_matrix/zeros_2/packed:output:0Ksequential_3/random_translation_1/translation_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������2>
<sequential_3/random_translation_1/translation_matrix/zeros_2�
@sequential_3/random_translation_1/translation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2B
@sequential_3/random_translation_1/translation_matrix/concat/axis�
;sequential_3/random_translation_1/translation_matrix/concatConcatV2Bsequential_3/random_translation_1/translation_matrix/ones:output:0Csequential_3/random_translation_1/translation_matrix/zeros:output:0<sequential_3/random_translation_1/translation_matrix/Neg:y:0Esequential_3/random_translation_1/translation_matrix/zeros_1:output:0Dsequential_3/random_translation_1/translation_matrix/ones_1:output:0>sequential_3/random_translation_1/translation_matrix/Neg_1:y:0Esequential_3/random_translation_1/translation_matrix/zeros_2:output:0Isequential_3/random_translation_1/translation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2=
;sequential_3/random_translation_1/translation_matrix/concat�
1sequential_3/random_translation_1/transform/ShapeShapeXsequential_3/random_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:23
1sequential_3/random_translation_1/transform/Shape�
?sequential_3/random_translation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2A
?sequential_3/random_translation_1/transform/strided_slice/stack�
Asequential_3/random_translation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2C
Asequential_3/random_translation_1/transform/strided_slice/stack_1�
Asequential_3/random_translation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2C
Asequential_3/random_translation_1/transform/strided_slice/stack_2�
9sequential_3/random_translation_1/transform/strided_sliceStridedSlice:sequential_3/random_translation_1/transform/Shape:output:0Hsequential_3/random_translation_1/transform/strided_slice/stack:output:0Jsequential_3/random_translation_1/transform/strided_slice/stack_1:output:0Jsequential_3/random_translation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2;
9sequential_3/random_translation_1/transform/strided_slice�
6sequential_3/random_translation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    28
6sequential_3/random_translation_1/transform/fill_value�
Fsequential_3/random_translation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Xsequential_3/random_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0Dsequential_3/random_translation_1/translation_matrix/concat:output:0Bsequential_3/random_translation_1/transform/strided_slice:output:0?sequential_3/random_translation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	NEAREST*
interpolation
BILINEAR2H
Fsequential_3/random_translation_1/transform/ImageProjectiveTransformV3�
 sequential_3/random_zoom_1/ShapeShape[sequential_3/random_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2"
 sequential_3/random_zoom_1/Shape�
.sequential_3/random_zoom_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 20
.sequential_3/random_zoom_1/strided_slice/stack�
0sequential_3/random_zoom_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:22
0sequential_3/random_zoom_1/strided_slice/stack_1�
0sequential_3/random_zoom_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:22
0sequential_3/random_zoom_1/strided_slice/stack_2�
(sequential_3/random_zoom_1/strided_sliceStridedSlice)sequential_3/random_zoom_1/Shape:output:07sequential_3/random_zoom_1/strided_slice/stack:output:09sequential_3/random_zoom_1/strided_slice/stack_1:output:09sequential_3/random_zoom_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2*
(sequential_3/random_zoom_1/strided_slice�
0sequential_3/random_zoom_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:22
0sequential_3/random_zoom_1/strided_slice_1/stack�
2sequential_3/random_zoom_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:24
2sequential_3/random_zoom_1/strided_slice_1/stack_1�
2sequential_3/random_zoom_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:24
2sequential_3/random_zoom_1/strided_slice_1/stack_2�
*sequential_3/random_zoom_1/strided_slice_1StridedSlice)sequential_3/random_zoom_1/Shape:output:09sequential_3/random_zoom_1/strided_slice_1/stack:output:0;sequential_3/random_zoom_1/strided_slice_1/stack_1:output:0;sequential_3/random_zoom_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2,
*sequential_3/random_zoom_1/strided_slice_1�
sequential_3/random_zoom_1/CastCast3sequential_3/random_zoom_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2!
sequential_3/random_zoom_1/Cast�
0sequential_3/random_zoom_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:22
0sequential_3/random_zoom_1/strided_slice_2/stack�
2sequential_3/random_zoom_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:24
2sequential_3/random_zoom_1/strided_slice_2/stack_1�
2sequential_3/random_zoom_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:24
2sequential_3/random_zoom_1/strided_slice_2/stack_2�
*sequential_3/random_zoom_1/strided_slice_2StridedSlice)sequential_3/random_zoom_1/Shape:output:09sequential_3/random_zoom_1/strided_slice_2/stack:output:0;sequential_3/random_zoom_1/strided_slice_2/stack_1:output:0;sequential_3/random_zoom_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2,
*sequential_3/random_zoom_1/strided_slice_2�
!sequential_3/random_zoom_1/Cast_1Cast3sequential_3/random_zoom_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2#
!sequential_3/random_zoom_1/Cast_1�
3sequential_3/random_zoom_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :25
3sequential_3/random_zoom_1/stateful_uniform/shape/1�
1sequential_3/random_zoom_1/stateful_uniform/shapePack1sequential_3/random_zoom_1/strided_slice:output:0<sequential_3/random_zoom_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:23
1sequential_3/random_zoom_1/stateful_uniform/shape�
/sequential_3/random_zoom_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *��L?21
/sequential_3/random_zoom_1/stateful_uniform/min�
/sequential_3/random_zoom_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���?21
/sequential_3/random_zoom_1/stateful_uniform/max�
Esequential_3/random_zoom_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2G
Esequential_3/random_zoom_1/stateful_uniform/StatefulUniform/algorithm�
;sequential_3/random_zoom_1/stateful_uniform/StatefulUniformStatefulUniformDsequential_3_random_zoom_1_stateful_uniform_statefuluniform_resourceNsequential_3/random_zoom_1/stateful_uniform/StatefulUniform/algorithm:output:0:sequential_3/random_zoom_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype02=
;sequential_3/random_zoom_1/stateful_uniform/StatefulUniform�
/sequential_3/random_zoom_1/stateful_uniform/subSub8sequential_3/random_zoom_1/stateful_uniform/max:output:08sequential_3/random_zoom_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 21
/sequential_3/random_zoom_1/stateful_uniform/sub�
/sequential_3/random_zoom_1/stateful_uniform/mulMulDsequential_3/random_zoom_1/stateful_uniform/StatefulUniform:output:03sequential_3/random_zoom_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������21
/sequential_3/random_zoom_1/stateful_uniform/mul�
+sequential_3/random_zoom_1/stateful_uniformAdd3sequential_3/random_zoom_1/stateful_uniform/mul:z:08sequential_3/random_zoom_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2-
+sequential_3/random_zoom_1/stateful_uniform�
&sequential_3/random_zoom_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2(
&sequential_3/random_zoom_1/concat/axis�
!sequential_3/random_zoom_1/concatConcatV2/sequential_3/random_zoom_1/stateful_uniform:z:0/sequential_3/random_zoom_1/stateful_uniform:z:0/sequential_3/random_zoom_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2#
!sequential_3/random_zoom_1/concat�
,sequential_3/random_zoom_1/zoom_matrix/ShapeShape*sequential_3/random_zoom_1/concat:output:0*
T0*
_output_shapes
:2.
,sequential_3/random_zoom_1/zoom_matrix/Shape�
:sequential_3/random_zoom_1/zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2<
:sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack�
<sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2>
<sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_1�
<sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2>
<sequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_2�
4sequential_3/random_zoom_1/zoom_matrix/strided_sliceStridedSlice5sequential_3/random_zoom_1/zoom_matrix/Shape:output:0Csequential_3/random_zoom_1/zoom_matrix/strided_slice/stack:output:0Esequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_1:output:0Esequential_3/random_zoom_1/zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask26
4sequential_3/random_zoom_1/zoom_matrix/strided_slice�
,sequential_3/random_zoom_1/zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2.
,sequential_3/random_zoom_1/zoom_matrix/sub/y�
*sequential_3/random_zoom_1/zoom_matrix/subSub%sequential_3/random_zoom_1/Cast_1:y:05sequential_3/random_zoom_1/zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 2,
*sequential_3/random_zoom_1/zoom_matrix/sub�
0sequential_3/random_zoom_1/zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @22
0sequential_3/random_zoom_1/zoom_matrix/truediv/y�
.sequential_3/random_zoom_1/zoom_matrix/truedivRealDiv.sequential_3/random_zoom_1/zoom_matrix/sub:z:09sequential_3/random_zoom_1/zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 20
.sequential_3/random_zoom_1/zoom_matrix/truediv�
<sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2>
<sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack�
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2@
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_1�
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2@
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_2�
6sequential_3/random_zoom_1/zoom_matrix/strided_slice_1StridedSlice*sequential_3/random_zoom_1/concat:output:0Esequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack:output:0Gsequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_1:output:0Gsequential_3/random_zoom_1/zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask28
6sequential_3/random_zoom_1/zoom_matrix/strided_slice_1�
.sequential_3/random_zoom_1/zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?20
.sequential_3/random_zoom_1/zoom_matrix/sub_1/x�
,sequential_3/random_zoom_1/zoom_matrix/sub_1Sub7sequential_3/random_zoom_1/zoom_matrix/sub_1/x:output:0?sequential_3/random_zoom_1/zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2.
,sequential_3/random_zoom_1/zoom_matrix/sub_1�
*sequential_3/random_zoom_1/zoom_matrix/mulMul2sequential_3/random_zoom_1/zoom_matrix/truediv:z:00sequential_3/random_zoom_1/zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:���������2,
*sequential_3/random_zoom_1/zoom_matrix/mul�
.sequential_3/random_zoom_1/zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?20
.sequential_3/random_zoom_1/zoom_matrix/sub_2/y�
,sequential_3/random_zoom_1/zoom_matrix/sub_2Sub#sequential_3/random_zoom_1/Cast:y:07sequential_3/random_zoom_1/zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2.
,sequential_3/random_zoom_1/zoom_matrix/sub_2�
2sequential_3/random_zoom_1/zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @24
2sequential_3/random_zoom_1/zoom_matrix/truediv_1/y�
0sequential_3/random_zoom_1/zoom_matrix/truediv_1RealDiv0sequential_3/random_zoom_1/zoom_matrix/sub_2:z:0;sequential_3/random_zoom_1/zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 22
0sequential_3/random_zoom_1/zoom_matrix/truediv_1�
<sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2>
<sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack�
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2@
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_1�
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2@
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_2�
6sequential_3/random_zoom_1/zoom_matrix/strided_slice_2StridedSlice*sequential_3/random_zoom_1/concat:output:0Esequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack:output:0Gsequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_1:output:0Gsequential_3/random_zoom_1/zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask28
6sequential_3/random_zoom_1/zoom_matrix/strided_slice_2�
.sequential_3/random_zoom_1/zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?20
.sequential_3/random_zoom_1/zoom_matrix/sub_3/x�
,sequential_3/random_zoom_1/zoom_matrix/sub_3Sub7sequential_3/random_zoom_1/zoom_matrix/sub_3/x:output:0?sequential_3/random_zoom_1/zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2.
,sequential_3/random_zoom_1/zoom_matrix/sub_3�
,sequential_3/random_zoom_1/zoom_matrix/mul_1Mul4sequential_3/random_zoom_1/zoom_matrix/truediv_1:z:00sequential_3/random_zoom_1/zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:���������2.
,sequential_3/random_zoom_1/zoom_matrix/mul_1�
<sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2>
<sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack�
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2@
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_1�
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2@
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_2�
6sequential_3/random_zoom_1/zoom_matrix/strided_slice_3StridedSlice*sequential_3/random_zoom_1/concat:output:0Esequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack:output:0Gsequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_1:output:0Gsequential_3/random_zoom_1/zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask28
6sequential_3/random_zoom_1/zoom_matrix/strided_slice_3�
2sequential_3/random_zoom_1/zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :24
2sequential_3/random_zoom_1/zoom_matrix/zeros/mul/y�
0sequential_3/random_zoom_1/zoom_matrix/zeros/mulMul=sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0;sequential_3/random_zoom_1/zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 22
0sequential_3/random_zoom_1/zoom_matrix/zeros/mul�
3sequential_3/random_zoom_1/zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�25
3sequential_3/random_zoom_1/zoom_matrix/zeros/Less/y�
1sequential_3/random_zoom_1/zoom_matrix/zeros/LessLess4sequential_3/random_zoom_1/zoom_matrix/zeros/mul:z:0<sequential_3/random_zoom_1/zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 23
1sequential_3/random_zoom_1/zoom_matrix/zeros/Less�
5sequential_3/random_zoom_1/zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :27
5sequential_3/random_zoom_1/zoom_matrix/zeros/packed/1�
3sequential_3/random_zoom_1/zoom_matrix/zeros/packedPack=sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0>sequential_3/random_zoom_1/zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:25
3sequential_3/random_zoom_1/zoom_matrix/zeros/packed�
2sequential_3/random_zoom_1/zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    24
2sequential_3/random_zoom_1/zoom_matrix/zeros/Const�
,sequential_3/random_zoom_1/zoom_matrix/zerosFill<sequential_3/random_zoom_1/zoom_matrix/zeros/packed:output:0;sequential_3/random_zoom_1/zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2.
,sequential_3/random_zoom_1/zoom_matrix/zeros�
4sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :26
4sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul/y�
2sequential_3/random_zoom_1/zoom_matrix/zeros_1/mulMul=sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0=sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 24
2sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul�
5sequential_3/random_zoom_1/zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�27
5sequential_3/random_zoom_1/zoom_matrix/zeros_1/Less/y�
3sequential_3/random_zoom_1/zoom_matrix/zeros_1/LessLess6sequential_3/random_zoom_1/zoom_matrix/zeros_1/mul:z:0>sequential_3/random_zoom_1/zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 25
3sequential_3/random_zoom_1/zoom_matrix/zeros_1/Less�
7sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :29
7sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed/1�
5sequential_3/random_zoom_1/zoom_matrix/zeros_1/packedPack=sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0@sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:27
5sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed�
4sequential_3/random_zoom_1/zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    26
4sequential_3/random_zoom_1/zoom_matrix/zeros_1/Const�
.sequential_3/random_zoom_1/zoom_matrix/zeros_1Fill>sequential_3/random_zoom_1/zoom_matrix/zeros_1/packed:output:0=sequential_3/random_zoom_1/zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������20
.sequential_3/random_zoom_1/zoom_matrix/zeros_1�
<sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2>
<sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack�
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2@
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_1�
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2@
>sequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_2�
6sequential_3/random_zoom_1/zoom_matrix/strided_slice_4StridedSlice*sequential_3/random_zoom_1/concat:output:0Esequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack:output:0Gsequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_1:output:0Gsequential_3/random_zoom_1/zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask28
6sequential_3/random_zoom_1/zoom_matrix/strided_slice_4�
4sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :26
4sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul/y�
2sequential_3/random_zoom_1/zoom_matrix/zeros_2/mulMul=sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0=sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 24
2sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul�
5sequential_3/random_zoom_1/zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�27
5sequential_3/random_zoom_1/zoom_matrix/zeros_2/Less/y�
3sequential_3/random_zoom_1/zoom_matrix/zeros_2/LessLess6sequential_3/random_zoom_1/zoom_matrix/zeros_2/mul:z:0>sequential_3/random_zoom_1/zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 25
3sequential_3/random_zoom_1/zoom_matrix/zeros_2/Less�
7sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :29
7sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed/1�
5sequential_3/random_zoom_1/zoom_matrix/zeros_2/packedPack=sequential_3/random_zoom_1/zoom_matrix/strided_slice:output:0@sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:27
5sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed�
4sequential_3/random_zoom_1/zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    26
4sequential_3/random_zoom_1/zoom_matrix/zeros_2/Const�
.sequential_3/random_zoom_1/zoom_matrix/zeros_2Fill>sequential_3/random_zoom_1/zoom_matrix/zeros_2/packed:output:0=sequential_3/random_zoom_1/zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������20
.sequential_3/random_zoom_1/zoom_matrix/zeros_2�
2sequential_3/random_zoom_1/zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :24
2sequential_3/random_zoom_1/zoom_matrix/concat/axis�
-sequential_3/random_zoom_1/zoom_matrix/concatConcatV2?sequential_3/random_zoom_1/zoom_matrix/strided_slice_3:output:05sequential_3/random_zoom_1/zoom_matrix/zeros:output:0.sequential_3/random_zoom_1/zoom_matrix/mul:z:07sequential_3/random_zoom_1/zoom_matrix/zeros_1:output:0?sequential_3/random_zoom_1/zoom_matrix/strided_slice_4:output:00sequential_3/random_zoom_1/zoom_matrix/mul_1:z:07sequential_3/random_zoom_1/zoom_matrix/zeros_2:output:0;sequential_3/random_zoom_1/zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2/
-sequential_3/random_zoom_1/zoom_matrix/concat�
*sequential_3/random_zoom_1/transform/ShapeShape[sequential_3/random_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2,
*sequential_3/random_zoom_1/transform/Shape�
8sequential_3/random_zoom_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2:
8sequential_3/random_zoom_1/transform/strided_slice/stack�
:sequential_3/random_zoom_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2<
:sequential_3/random_zoom_1/transform/strided_slice/stack_1�
:sequential_3/random_zoom_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2<
:sequential_3/random_zoom_1/transform/strided_slice/stack_2�
2sequential_3/random_zoom_1/transform/strided_sliceStridedSlice3sequential_3/random_zoom_1/transform/Shape:output:0Asequential_3/random_zoom_1/transform/strided_slice/stack:output:0Csequential_3/random_zoom_1/transform/strided_slice/stack_1:output:0Csequential_3/random_zoom_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:24
2sequential_3/random_zoom_1/transform/strided_slice�
/sequential_3/random_zoom_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    21
/sequential_3/random_zoom_1/transform/fill_value�
?sequential_3/random_zoom_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3[sequential_3/random_translation_1/transform/ImageProjectiveTransformV3:transformed_images:06sequential_3/random_zoom_1/zoom_matrix/concat:output:0;sequential_3/random_zoom_1/transform/strided_slice:output:08sequential_3/random_zoom_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2A
?sequential_3/random_zoom_1/transform/ImageProjectiveTransformV3m
rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2
rescaling_3/Cast/xq
rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rescaling_3/Cast_1/x�
rescaling_3/mulMulTsequential_3/random_zoom_1/transform/ImageProjectiveTransformV3:transformed_images:0rescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/mul�
rescaling_3/addAddV2rescaling_3/mul:z:0rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/add�
conv2d_3/Conv2D/ReadVariableOpReadVariableOp'conv2d_3_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02 
conv2d_3/Conv2D/ReadVariableOp�
conv2d_3/Conv2DConv2Drescaling_3/add:z:0&conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&*
paddingSAME*
strides
2
conv2d_3/Conv2D�
conv2d_3/BiasAdd/ReadVariableOpReadVariableOp(conv2d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv2d_3/BiasAdd/ReadVariableOp�
conv2d_3/BiasAddBiasAddconv2d_3/Conv2D:output:0'conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&2
conv2d_3/BiasAdd{
conv2d_3/ReluReluconv2d_3/BiasAdd:output:0*
T0*/
_output_shapes
:���������&2
conv2d_3/Reluw
dropout_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout_4/dropout/Const�
dropout_4/dropout/MulMulconv2d_3/Relu:activations:0 dropout_4/dropout/Const:output:0*
T0*/
_output_shapes
:���������&2
dropout_4/dropout/Mul}
dropout_4/dropout/ShapeShapeconv2d_3/Relu:activations:0*
T0*
_output_shapes
:2
dropout_4/dropout/Shape�
.dropout_4/dropout/random_uniform/RandomUniformRandomUniform dropout_4/dropout/Shape:output:0*
T0*/
_output_shapes
:���������&*
dtype020
.dropout_4/dropout/random_uniform/RandomUniform�
 dropout_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2"
 dropout_4/dropout/GreaterEqual/y�
dropout_4/dropout/GreaterEqualGreaterEqual7dropout_4/dropout/random_uniform/RandomUniform:output:0)dropout_4/dropout/GreaterEqual/y:output:0*
T0*/
_output_shapes
:���������&2 
dropout_4/dropout/GreaterEqual�
dropout_4/dropout/CastCast"dropout_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*/
_output_shapes
:���������&2
dropout_4/dropout/Cast�
dropout_4/dropout/Mul_1Muldropout_4/dropout/Mul:z:0dropout_4/dropout/Cast:y:0*
T0*/
_output_shapes
:���������&2
dropout_4/dropout/Mul_1�
max_pooling2d_3/MaxPoolMaxPooldropout_4/dropout/Mul_1:z:0*/
_output_shapes
:���������*
ksize
*
paddingVALID*
strides
2
max_pooling2d_3/MaxPool�
conv2d_4/Conv2D/ReadVariableOpReadVariableOp'conv2d_4_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02 
conv2d_4/Conv2D/ReadVariableOp�
conv2d_4/Conv2DConv2D max_pooling2d_3/MaxPool:output:0&conv2d_4/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������*
paddingSAME*
strides
2
conv2d_4/Conv2D�
conv2d_4/BiasAdd/ReadVariableOpReadVariableOp(conv2d_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv2d_4/BiasAdd/ReadVariableOp�
conv2d_4/BiasAddBiasAddconv2d_4/Conv2D:output:0'conv2d_4/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������2
conv2d_4/BiasAdd{
conv2d_4/ReluReluconv2d_4/BiasAdd:output:0*
T0*/
_output_shapes
:���������2
conv2d_4/Reluw
dropout_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout_5/dropout/Const�
dropout_5/dropout/MulMulconv2d_4/Relu:activations:0 dropout_5/dropout/Const:output:0*
T0*/
_output_shapes
:���������2
dropout_5/dropout/Mul}
dropout_5/dropout/ShapeShapeconv2d_4/Relu:activations:0*
T0*
_output_shapes
:2
dropout_5/dropout/Shape�
.dropout_5/dropout/random_uniform/RandomUniformRandomUniform dropout_5/dropout/Shape:output:0*
T0*/
_output_shapes
:���������*
dtype020
.dropout_5/dropout/random_uniform/RandomUniform�
 dropout_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2"
 dropout_5/dropout/GreaterEqual/y�
dropout_5/dropout/GreaterEqualGreaterEqual7dropout_5/dropout/random_uniform/RandomUniform:output:0)dropout_5/dropout/GreaterEqual/y:output:0*
T0*/
_output_shapes
:���������2 
dropout_5/dropout/GreaterEqual�
dropout_5/dropout/CastCast"dropout_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*/
_output_shapes
:���������2
dropout_5/dropout/Cast�
dropout_5/dropout/Mul_1Muldropout_5/dropout/Mul:z:0dropout_5/dropout/Cast:y:0*
T0*/
_output_shapes
:���������2
dropout_5/dropout/Mul_1�
max_pooling2d_4/MaxPoolMaxPooldropout_5/dropout/Mul_1:z:0*/
_output_shapes
:���������	*
ksize
*
paddingVALID*
strides
2
max_pooling2d_4/MaxPool�
conv2d_5/Conv2D/ReadVariableOpReadVariableOp'conv2d_5_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype02 
conv2d_5/Conv2D/ReadVariableOp�
conv2d_5/Conv2DConv2D max_pooling2d_4/MaxPool:output:0&conv2d_5/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 *
paddingSAME*
strides
2
conv2d_5/Conv2D�
conv2d_5/BiasAdd/ReadVariableOpReadVariableOp(conv2d_5_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02!
conv2d_5/BiasAdd/ReadVariableOp�
conv2d_5/BiasAddBiasAddconv2d_5/Conv2D:output:0'conv2d_5/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 2
conv2d_5/BiasAdd{
conv2d_5/ReluReluconv2d_5/BiasAdd:output:0*
T0*/
_output_shapes
:���������	 2
conv2d_5/Relu�
max_pooling2d_5/MaxPoolMaxPoolconv2d_5/Relu:activations:0*/
_output_shapes
:��������� *
ksize
*
paddingVALID*
strides
2
max_pooling2d_5/MaxPools
flatten_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"�����  2
flatten_1/Const�
flatten_1/ReshapeReshape max_pooling2d_5/MaxPool:output:0flatten_1/Const:output:0*
T0*(
_output_shapes
:����������2
flatten_1/Reshape�
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype02
dense_3/MatMul/ReadVariableOp�
dense_3/MatMulMatMulflatten_1/Reshape:output:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
dense_3/MatMul�
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype02 
dense_3/BiasAdd/ReadVariableOp�
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
dense_3/BiasAddq
dense_3/ReluReludense_3/BiasAdd:output:0*
T0*(
_output_shapes
:����������2
dense_3/Reluw
dropout_6/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_6/dropout/Const�
dropout_6/dropout/MulMuldense_3/Relu:activations:0 dropout_6/dropout/Const:output:0*
T0*(
_output_shapes
:����������2
dropout_6/dropout/Mul|
dropout_6/dropout/ShapeShapedense_3/Relu:activations:0*
T0*
_output_shapes
:2
dropout_6/dropout/Shape�
.dropout_6/dropout/random_uniform/RandomUniformRandomUniform dropout_6/dropout/Shape:output:0*
T0*(
_output_shapes
:����������*
dtype020
.dropout_6/dropout/random_uniform/RandomUniform�
 dropout_6/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2"
 dropout_6/dropout/GreaterEqual/y�
dropout_6/dropout/GreaterEqualGreaterEqual7dropout_6/dropout/random_uniform/RandomUniform:output:0)dropout_6/dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:����������2 
dropout_6/dropout/GreaterEqual�
dropout_6/dropout/CastCast"dropout_6/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:����������2
dropout_6/dropout/Cast�
dropout_6/dropout/Mul_1Muldropout_6/dropout/Mul:z:0dropout_6/dropout/Cast:y:0*
T0*(
_output_shapes
:����������2
dropout_6/dropout/Mul_1�
dense_4/MatMul/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype02
dense_4/MatMul/ReadVariableOp�
dense_4/MatMulMatMuldropout_6/dropout/Mul_1:z:0%dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_4/MatMul�
dense_4/BiasAdd/ReadVariableOpReadVariableOp'dense_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_4/BiasAdd/ReadVariableOp�
dense_4/BiasAddBiasAdddense_4/MatMul:product:0&dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_4/BiasAddp
dense_4/ReluReludense_4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
dense_4/Reluw
dropout_7/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout_7/dropout/Const�
dropout_7/dropout/MulMuldense_4/Relu:activations:0 dropout_7/dropout/Const:output:0*
T0*'
_output_shapes
:���������2
dropout_7/dropout/Mul|
dropout_7/dropout/ShapeShapedense_4/Relu:activations:0*
T0*
_output_shapes
:2
dropout_7/dropout/Shape�
.dropout_7/dropout/random_uniform/RandomUniformRandomUniform dropout_7/dropout/Shape:output:0*
T0*'
_output_shapes
:���������*
dtype020
.dropout_7/dropout/random_uniform/RandomUniform�
 dropout_7/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���>2"
 dropout_7/dropout/GreaterEqual/y�
dropout_7/dropout/GreaterEqualGreaterEqual7dropout_7/dropout/random_uniform/RandomUniform:output:0)dropout_7/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:���������2 
dropout_7/dropout/GreaterEqual�
dropout_7/dropout/CastCast"dropout_7/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:���������2
dropout_7/dropout/Cast�
dropout_7/dropout/Mul_1Muldropout_7/dropout/Mul:z:0dropout_7/dropout/Cast:y:0*
T0*'
_output_shapes
:���������2
dropout_7/dropout/Mul_1�
dense_5/MatMul/ReadVariableOpReadVariableOp&dense_5_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
dense_5/MatMul/ReadVariableOp�
dense_5/MatMulMatMuldropout_7/dropout/Mul_1:z:0%dense_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_5/MatMul�
dense_5/BiasAdd/ReadVariableOpReadVariableOp'dense_5_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_5/BiasAdd/ReadVariableOp�
dense_5/BiasAddBiasAdddense_5/MatMul:product:0&dense_5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_5/BiasAddy
dense_5/SigmoidSigmoiddense_5/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
dense_5/Sigmoid�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentitydense_5/Sigmoid:y:0 ^conv2d_3/BiasAdd/ReadVariableOp^conv2d_3/Conv2D/ReadVariableOp ^conv2d_4/BiasAdd/ReadVariableOp^conv2d_4/Conv2D/ReadVariableOp ^conv2d_5/BiasAdd/ReadVariableOp^conv2d_5/Conv2D/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp^dense_4/BiasAdd/ReadVariableOp^dense_4/MatMul/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp^dense_5/BiasAdd/ReadVariableOp^dense_5/MatMul/ReadVariableOp@^sequential_3/random_rotation_1/stateful_uniform/StatefulUniformC^sequential_3/random_translation_1/stateful_uniform/StatefulUniformE^sequential_3/random_translation_1/stateful_uniform_1/StatefulUniform<^sequential_3/random_zoom_1/stateful_uniform/StatefulUniform*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::2B
conv2d_3/BiasAdd/ReadVariableOpconv2d_3/BiasAdd/ReadVariableOp2@
conv2d_3/Conv2D/ReadVariableOpconv2d_3/Conv2D/ReadVariableOp2B
conv2d_4/BiasAdd/ReadVariableOpconv2d_4/BiasAdd/ReadVariableOp2@
conv2d_4/Conv2D/ReadVariableOpconv2d_4/Conv2D/ReadVariableOp2B
conv2d_5/BiasAdd/ReadVariableOpconv2d_5/BiasAdd/ReadVariableOp2@
conv2d_5/Conv2D/ReadVariableOpconv2d_5/Conv2D/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2@
dense_4/BiasAdd/ReadVariableOpdense_4/BiasAdd/ReadVariableOp2>
dense_4/MatMul/ReadVariableOpdense_4/MatMul/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2@
dense_5/BiasAdd/ReadVariableOpdense_5/BiasAdd/ReadVariableOp2>
dense_5/MatMul/ReadVariableOpdense_5/MatMul/ReadVariableOp2�
?sequential_3/random_rotation_1/stateful_uniform/StatefulUniform?sequential_3/random_rotation_1/stateful_uniform/StatefulUniform2�
Bsequential_3/random_translation_1/stateful_uniform/StatefulUniformBsequential_3/random_translation_1/stateful_uniform/StatefulUniform2�
Dsequential_3/random_translation_1/stateful_uniform_1/StatefulUniformDsequential_3/random_translation_1/stateful_uniform_1/StatefulUniform2z
;sequential_3/random_zoom_1/stateful_uniform/StatefulUniform;sequential_3/random_zoom_1/stateful_uniform/StatefulUniform:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
f
J__inference_max_pooling2d_5_layer_call_and_return_conditional_losses_31341

inputs
identity�
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4������������������������������������*
ksize
*
paddingVALID*
strides
2	
MaxPool�
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4������������������������������������2

Identity"
identityIdentity:output:0*I
_input_shapes8
6:4������������������������������������:r n
J
_output_shapes8
6:4������������������������������������
 
_user_specified_nameinputs
�	
�
,__inference_sequential_4_layer_call_fn_33453

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_4_layer_call_and_return_conditional_losses_319232
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�

�
,__inference_sequential_5_layer_call_fn_32846

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_5_layer_call_and_return_conditional_losses_321682
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
�
__inference_loss_fn_1_34096=
9dense_4_kernel_regularizer_square_readvariableop_resource
identity��0dense_4/kernel/Regularizer/Square/ReadVariableOp�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_4_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity"dense_4/kernel/Regularizer/mul:z:01^dense_4/kernel/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*
_input_shapes
:2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp
�	
�
,__inference_sequential_4_layer_call_fn_31950
sequential_3_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallsequential_3_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_4_layer_call_and_return_conditional_losses_319232
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_3_input
�#
�
G__inference_sequential_5_layer_call_and_return_conditional_losses_32075
sequential_4_input
sequential_4_32018
sequential_4_32020
sequential_4_32022
sequential_4_32024
sequential_4_32026
sequential_4_32028
sequential_4_32030
sequential_4_32032
sequential_4_32034
sequential_4_32036
sequential_4_32038
sequential_4_32040
sequential_4_32042
sequential_4_32044
sequential_4_32046
identity��0dense_3/kernel/Regularizer/Square/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�$sequential_4/StatefulPartitionedCall�
$sequential_4/StatefulPartitionedCallStatefulPartitionedCallsequential_4_inputsequential_4_32018sequential_4_32020sequential_4_32022sequential_4_32024sequential_4_32026sequential_4_32028sequential_4_32030sequential_4_32032sequential_4_32034sequential_4_32036sequential_4_32038sequential_4_32040sequential_4_32042sequential_4_32044sequential_4_32046*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_4_layer_call_and_return_conditional_losses_318292&
$sequential_4/StatefulPartitionedCall�
softmax_1/PartitionedCallPartitionedCall-sequential_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_softmax_1_layer_call_and_return_conditional_losses_320542
softmax_1/PartitionedCall�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_32036* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_32040*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity"softmax_1/PartitionedCall:output:01^dense_3/kernel/Regularizer/Square/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp%^sequential_4/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_4_input
�
�
,__inference_sequential_3_layer_call_fn_31302
random_flip_1_input
unknown
	unknown_0
	unknown_1
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallrandom_flip_1_inputunknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_312932
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*:
_input_shapes)
':���������&:::22
StatefulPartitionedCallStatefulPartitionedCall:d `
/
_output_shapes
:���������&
-
_user_specified_namerandom_flip_1_input
�
�
B__inference_dense_4_layer_call_and_return_conditional_losses_34018

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������2
Relu�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�#
�
G__inference_sequential_5_layer_call_and_return_conditional_losses_32168

inputs
sequential_4_32123
sequential_4_32125
sequential_4_32127
sequential_4_32129
sequential_4_32131
sequential_4_32133
sequential_4_32135
sequential_4_32137
sequential_4_32139
sequential_4_32141
sequential_4_32143
sequential_4_32145
sequential_4_32147
sequential_4_32149
sequential_4_32151
identity��0dense_3/kernel/Regularizer/Square/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�$sequential_4/StatefulPartitionedCall�
$sequential_4/StatefulPartitionedCallStatefulPartitionedCallinputssequential_4_32123sequential_4_32125sequential_4_32127sequential_4_32129sequential_4_32131sequential_4_32133sequential_4_32135sequential_4_32137sequential_4_32139sequential_4_32141sequential_4_32143sequential_4_32145sequential_4_32147sequential_4_32149sequential_4_32151*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_4_layer_call_and_return_conditional_losses_318292&
$sequential_4/StatefulPartitionedCall�
softmax_1/PartitionedCallPartitionedCall-sequential_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_softmax_1_layer_call_and_return_conditional_losses_320542
softmax_1/PartitionedCall�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_32141* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_4_32145*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity"softmax_1/PartitionedCall:output:01^dense_3/kernel/Regularizer/Square/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp%^sequential_4/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2L
$sequential_4/StatefulPartitionedCall$sequential_4/StatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
E
)__inference_dropout_6_layer_call_fn_33995

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_6_layer_call_and_return_conditional_losses_315862
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*'
_input_shapes
:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�S
�
G__inference_sequential_4_layer_call_and_return_conditional_losses_31923

inputs
conv2d_3_31872
conv2d_3_31874
conv2d_4_31879
conv2d_4_31881
conv2d_5_31886
conv2d_5_31888
dense_3_31893
dense_3_31895
dense_4_31899
dense_4_31901
dense_5_31905
dense_5_31907
identity�� conv2d_3/StatefulPartitionedCall� conv2d_4/StatefulPartitionedCall� conv2d_5/StatefulPartitionedCall�dense_3/StatefulPartitionedCall�0dense_3/kernel/Regularizer/Square/ReadVariableOp�dense_4/StatefulPartitionedCall�0dense_4/kernel/Regularizer/Square/ReadVariableOp�dense_5/StatefulPartitionedCall�
sequential_3/PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_313082
sequential_3/PartitionedCallm
rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2
rescaling_3/Cast/xq
rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rescaling_3/Cast_1/x�
rescaling_3/mulMul%sequential_3/PartitionedCall:output:0rescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/mul�
rescaling_3/addAddV2rescaling_3/mul:z:0rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/add�
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCallrescaling_3/add:z:0conv2d_3_31872conv2d_3_31874*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_313892"
 conv2d_3/StatefulPartitionedCall�
dropout_4/PartitionedCallPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_4_layer_call_and_return_conditional_losses_314222
dropout_4/PartitionedCall�
max_pooling2d_3/PartitionedCallPartitionedCall"dropout_4/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_313172!
max_pooling2d_3/PartitionedCall�
 conv2d_4/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_3/PartitionedCall:output:0conv2d_4_31879conv2d_4_31881*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_314472"
 conv2d_4/StatefulPartitionedCall�
dropout_5/PartitionedCallPartitionedCall)conv2d_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_5_layer_call_and_return_conditional_losses_314802
dropout_5/PartitionedCall�
max_pooling2d_4/PartitionedCallPartitionedCall"dropout_5/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_4_layer_call_and_return_conditional_losses_313292!
max_pooling2d_4/PartitionedCall�
 conv2d_5/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_4/PartitionedCall:output:0conv2d_5_31886conv2d_5_31888*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_5_layer_call_and_return_conditional_losses_315052"
 conv2d_5/StatefulPartitionedCall�
max_pooling2d_5/PartitionedCallPartitionedCall)conv2d_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_5_layer_call_and_return_conditional_losses_313412!
max_pooling2d_5/PartitionedCall�
flatten_1/PartitionedCallPartitionedCall(max_pooling2d_5/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_flatten_1_layer_call_and_return_conditional_losses_315282
flatten_1/PartitionedCall�
dense_3/StatefulPartitionedCallStatefulPartitionedCall"flatten_1/PartitionedCall:output:0dense_3_31893dense_3_31895*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_315532!
dense_3/StatefulPartitionedCall�
dropout_6/PartitionedCallPartitionedCall(dense_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_6_layer_call_and_return_conditional_losses_315862
dropout_6/PartitionedCall�
dense_4/StatefulPartitionedCallStatefulPartitionedCall"dropout_6/PartitionedCall:output:0dense_4_31899dense_4_31901*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_4_layer_call_and_return_conditional_losses_316162!
dense_4/StatefulPartitionedCall�
dropout_7/PartitionedCallPartitionedCall(dense_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_7_layer_call_and_return_conditional_losses_316492
dropout_7/PartitionedCall�
dense_5/StatefulPartitionedCallStatefulPartitionedCall"dropout_7/PartitionedCall:output:0dense_5_31905dense_5_31907*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_5_layer_call_and_return_conditional_losses_316732!
dense_5/StatefulPartitionedCall�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_31893* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_4_31899*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity(dense_5/StatefulPartitionedCall:output:0!^conv2d_3/StatefulPartitionedCall!^conv2d_4/StatefulPartitionedCall!^conv2d_5/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall1^dense_3/kernel/Regularizer/Square/ReadVariableOp ^dense_4/StatefulPartitionedCall1^dense_4/kernel/Regularizer/Square/ReadVariableOp ^dense_5/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2D
 conv2d_4/StatefulPartitionedCall conv2d_4/StatefulPartitionedCall2D
 conv2d_5/StatefulPartitionedCall conv2d_5/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2B
dense_5/StatefulPartitionedCalldense_5/StatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
c
D__inference_dropout_6_layer_call_and_return_conditional_losses_33980

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:����������2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:����������*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:����������2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:����������2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:����������2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*'
_input_shapes
:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�	
�
#__inference_signature_wrapper_32315
sequential_4_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallsequential_4_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *)
f$R"
 __inference__wrapped_model_306292
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_4_input
�
c
G__inference_sequential_3_layer_call_and_return_conditional_losses_31308

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
c
D__inference_dropout_6_layer_call_and_return_conditional_losses_31581

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:����������2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:����������*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:����������2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:����������2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:����������2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*'
_input_shapes
:����������:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
,__inference_sequential_3_layer_call_fn_33806

inputs
unknown
	unknown_0
	unknown_1
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_312932
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*:
_input_shapes)
':���������&:::22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
`
D__inference_flatten_1_layer_call_and_return_conditional_losses_33931

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"�����  2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:����������2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*.
_input_shapes
:��������� :W S
/
_output_shapes
:��������� 
 
_user_specified_nameinputs
׾
�
G__inference_sequential_3_layer_call_and_return_conditional_losses_30958
random_flip_1_input?
;random_rotation_1_stateful_uniform_statefuluniform_resourceB
>random_translation_1_stateful_uniform_statefuluniform_resource;
7random_zoom_1_stateful_uniform_statefuluniform_resource
identity��2random_rotation_1/stateful_uniform/StatefulUniform�5random_translation_1/stateful_uniform/StatefulUniform�7random_translation_1/stateful_uniform_1/StatefulUniform�.random_zoom_1/stateful_uniform/StatefulUniform�
7random_flip_1/random_flip_left_right/control_dependencyIdentityrandom_flip_1_input*
T0*&
_class
loc:@random_flip_1_input*/
_output_shapes
:���������&29
7random_flip_1/random_flip_left_right/control_dependency�
*random_flip_1/random_flip_left_right/ShapeShape@random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:2,
*random_flip_1/random_flip_left_right/Shape�
8random_flip_1/random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2:
8random_flip_1/random_flip_left_right/strided_slice/stack�
:random_flip_1/random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2<
:random_flip_1/random_flip_left_right/strided_slice/stack_1�
:random_flip_1/random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2<
:random_flip_1/random_flip_left_right/strided_slice/stack_2�
2random_flip_1/random_flip_left_right/strided_sliceStridedSlice3random_flip_1/random_flip_left_right/Shape:output:0Arandom_flip_1/random_flip_left_right/strided_slice/stack:output:0Crandom_flip_1/random_flip_left_right/strided_slice/stack_1:output:0Crandom_flip_1/random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask24
2random_flip_1/random_flip_left_right/strided_slice�
9random_flip_1/random_flip_left_right/random_uniform/shapePack;random_flip_1/random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2;
9random_flip_1/random_flip_left_right/random_uniform/shape�
7random_flip_1/random_flip_left_right/random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    29
7random_flip_1/random_flip_left_right/random_uniform/min�
7random_flip_1/random_flip_left_right/random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  �?29
7random_flip_1/random_flip_left_right/random_uniform/max�
Arandom_flip_1/random_flip_left_right/random_uniform/RandomUniformRandomUniformBrandom_flip_1/random_flip_left_right/random_uniform/shape:output:0*
T0*#
_output_shapes
:���������*
dtype02C
Arandom_flip_1/random_flip_left_right/random_uniform/RandomUniform�
7random_flip_1/random_flip_left_right/random_uniform/MulMulJrandom_flip_1/random_flip_left_right/random_uniform/RandomUniform:output:0@random_flip_1/random_flip_left_right/random_uniform/max:output:0*
T0*#
_output_shapes
:���������29
7random_flip_1/random_flip_left_right/random_uniform/Mul�
4random_flip_1/random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/1�
4random_flip_1/random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/2�
4random_flip_1/random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/3�
2random_flip_1/random_flip_left_right/Reshape/shapePack;random_flip_1/random_flip_left_right/strided_slice:output:0=random_flip_1/random_flip_left_right/Reshape/shape/1:output:0=random_flip_1/random_flip_left_right/Reshape/shape/2:output:0=random_flip_1/random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:24
2random_flip_1/random_flip_left_right/Reshape/shape�
,random_flip_1/random_flip_left_right/ReshapeReshape;random_flip_1/random_flip_left_right/random_uniform/Mul:z:0;random_flip_1/random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:���������2.
,random_flip_1/random_flip_left_right/Reshape�
*random_flip_1/random_flip_left_right/RoundRound5random_flip_1/random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:���������2,
*random_flip_1/random_flip_left_right/Round�
3random_flip_1/random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:25
3random_flip_1/random_flip_left_right/ReverseV2/axis�
.random_flip_1/random_flip_left_right/ReverseV2	ReverseV2@random_flip_1/random_flip_left_right/control_dependency:output:0<random_flip_1/random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:���������&20
.random_flip_1/random_flip_left_right/ReverseV2�
(random_flip_1/random_flip_left_right/mulMul.random_flip_1/random_flip_left_right/Round:y:07random_flip_1/random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:���������&2*
(random_flip_1/random_flip_left_right/mul�
*random_flip_1/random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2,
*random_flip_1/random_flip_left_right/sub/x�
(random_flip_1/random_flip_left_right/subSub3random_flip_1/random_flip_left_right/sub/x:output:0.random_flip_1/random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:���������2*
(random_flip_1/random_flip_left_right/sub�
*random_flip_1/random_flip_left_right/mul_1Mul,random_flip_1/random_flip_left_right/sub:z:0@random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:���������&2,
*random_flip_1/random_flip_left_right/mul_1�
(random_flip_1/random_flip_left_right/addAddV2,random_flip_1/random_flip_left_right/mul:z:0.random_flip_1/random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:���������&2*
(random_flip_1/random_flip_left_right/add�
random_rotation_1/ShapeShape,random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:2
random_rotation_1/Shape�
%random_rotation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2'
%random_rotation_1/strided_slice/stack�
'random_rotation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice/stack_1�
'random_rotation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice/stack_2�
random_rotation_1/strided_sliceStridedSlice random_rotation_1/Shape:output:0.random_rotation_1/strided_slice/stack:output:00random_rotation_1/strided_slice/stack_1:output:00random_rotation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2!
random_rotation_1/strided_slice�
'random_rotation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice_1/stack�
)random_rotation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_1/stack_1�
)random_rotation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_1/stack_2�
!random_rotation_1/strided_slice_1StridedSlice random_rotation_1/Shape:output:00random_rotation_1/strided_slice_1/stack:output:02random_rotation_1/strided_slice_1/stack_1:output:02random_rotation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2#
!random_rotation_1/strided_slice_1�
random_rotation_1/CastCast*random_rotation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_rotation_1/Cast�
'random_rotation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice_2/stack�
)random_rotation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_2/stack_1�
)random_rotation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_2/stack_2�
!random_rotation_1/strided_slice_2StridedSlice random_rotation_1/Shape:output:00random_rotation_1/strided_slice_2/stack:output:02random_rotation_1/strided_slice_2/stack_1:output:02random_rotation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2#
!random_rotation_1/strided_slice_2�
random_rotation_1/Cast_1Cast*random_rotation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_rotation_1/Cast_1�
(random_rotation_1/stateful_uniform/shapePack(random_rotation_1/strided_slice:output:0*
N*
T0*
_output_shapes
:2*
(random_rotation_1/stateful_uniform/shape�
&random_rotation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *   �2(
&random_rotation_1/stateful_uniform/min�
&random_rotation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2(
&random_rotation_1/stateful_uniform/max�
<random_rotation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2>
<random_rotation_1/stateful_uniform/StatefulUniform/algorithm�
2random_rotation_1/stateful_uniform/StatefulUniformStatefulUniform;random_rotation_1_stateful_uniform_statefuluniform_resourceErandom_rotation_1/stateful_uniform/StatefulUniform/algorithm:output:01random_rotation_1/stateful_uniform/shape:output:0*#
_output_shapes
:���������*
shape_dtype024
2random_rotation_1/stateful_uniform/StatefulUniform�
&random_rotation_1/stateful_uniform/subSub/random_rotation_1/stateful_uniform/max:output:0/random_rotation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2(
&random_rotation_1/stateful_uniform/sub�
&random_rotation_1/stateful_uniform/mulMul;random_rotation_1/stateful_uniform/StatefulUniform:output:0*random_rotation_1/stateful_uniform/sub:z:0*
T0*#
_output_shapes
:���������2(
&random_rotation_1/stateful_uniform/mul�
"random_rotation_1/stateful_uniformAdd*random_rotation_1/stateful_uniform/mul:z:0/random_rotation_1/stateful_uniform/min:output:0*
T0*#
_output_shapes
:���������2$
"random_rotation_1/stateful_uniform�
'random_rotation_1/rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2)
'random_rotation_1/rotation_matrix/sub/y�
%random_rotation_1/rotation_matrix/subSubrandom_rotation_1/Cast_1:y:00random_rotation_1/rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 2'
%random_rotation_1/rotation_matrix/sub�
%random_rotation_1/rotation_matrix/CosCos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Cos�
)random_rotation_1/rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_1/y�
'random_rotation_1/rotation_matrix/sub_1Subrandom_rotation_1/Cast_1:y:02random_rotation_1/rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_1�
%random_rotation_1/rotation_matrix/mulMul)random_rotation_1/rotation_matrix/Cos:y:0+random_rotation_1/rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/mul�
%random_rotation_1/rotation_matrix/SinSin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Sin�
)random_rotation_1/rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_2/y�
'random_rotation_1/rotation_matrix/sub_2Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_2�
'random_rotation_1/rotation_matrix/mul_1Mul)random_rotation_1/rotation_matrix/Sin:y:0+random_rotation_1/rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_1�
'random_rotation_1/rotation_matrix/sub_3Sub)random_rotation_1/rotation_matrix/mul:z:0+random_rotation_1/rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_3�
'random_rotation_1/rotation_matrix/sub_4Sub)random_rotation_1/rotation_matrix/sub:z:0+random_rotation_1/rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_4�
+random_rotation_1/rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2-
+random_rotation_1/rotation_matrix/truediv/y�
)random_rotation_1/rotation_matrix/truedivRealDiv+random_rotation_1/rotation_matrix/sub_4:z:04random_rotation_1/rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:���������2+
)random_rotation_1/rotation_matrix/truediv�
)random_rotation_1/rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_5/y�
'random_rotation_1/rotation_matrix/sub_5Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_5�
'random_rotation_1/rotation_matrix/Sin_1Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_1�
)random_rotation_1/rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_6/y�
'random_rotation_1/rotation_matrix/sub_6Subrandom_rotation_1/Cast_1:y:02random_rotation_1/rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_6�
'random_rotation_1/rotation_matrix/mul_2Mul+random_rotation_1/rotation_matrix/Sin_1:y:0+random_rotation_1/rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_2�
'random_rotation_1/rotation_matrix/Cos_1Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_1�
)random_rotation_1/rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_7/y�
'random_rotation_1/rotation_matrix/sub_7Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_7�
'random_rotation_1/rotation_matrix/mul_3Mul+random_rotation_1/rotation_matrix/Cos_1:y:0+random_rotation_1/rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_3�
%random_rotation_1/rotation_matrix/addAddV2+random_rotation_1/rotation_matrix/mul_2:z:0+random_rotation_1/rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/add�
'random_rotation_1/rotation_matrix/sub_8Sub+random_rotation_1/rotation_matrix/sub_5:z:0)random_rotation_1/rotation_matrix/add:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_8�
-random_rotation_1/rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2/
-random_rotation_1/rotation_matrix/truediv_1/y�
+random_rotation_1/rotation_matrix/truediv_1RealDiv+random_rotation_1/rotation_matrix/sub_8:z:06random_rotation_1/rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:���������2-
+random_rotation_1/rotation_matrix/truediv_1�
'random_rotation_1/rotation_matrix/ShapeShape&random_rotation_1/stateful_uniform:z:0*
T0*
_output_shapes
:2)
'random_rotation_1/rotation_matrix/Shape�
5random_rotation_1/rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 27
5random_rotation_1/rotation_matrix/strided_slice/stack�
7random_rotation_1/rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:29
7random_rotation_1/rotation_matrix/strided_slice/stack_1�
7random_rotation_1/rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:29
7random_rotation_1/rotation_matrix/strided_slice/stack_2�
/random_rotation_1/rotation_matrix/strided_sliceStridedSlice0random_rotation_1/rotation_matrix/Shape:output:0>random_rotation_1/rotation_matrix/strided_slice/stack:output:0@random_rotation_1/rotation_matrix/strided_slice/stack_1:output:0@random_rotation_1/rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask21
/random_rotation_1/rotation_matrix/strided_slice�
'random_rotation_1/rotation_matrix/Cos_2Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_2�
7random_rotation_1/rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_1/stack�
9random_rotation_1/rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_1/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_1/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_1StridedSlice+random_rotation_1/rotation_matrix/Cos_2:y:0@random_rotation_1/rotation_matrix/strided_slice_1/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_1/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_1�
'random_rotation_1/rotation_matrix/Sin_2Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_2�
7random_rotation_1/rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_2/stack�
9random_rotation_1/rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_2/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_2/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_2StridedSlice+random_rotation_1/rotation_matrix/Sin_2:y:0@random_rotation_1/rotation_matrix/strided_slice_2/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_2/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_2�
%random_rotation_1/rotation_matrix/NegNeg:random_rotation_1/rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Neg�
7random_rotation_1/rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_3/stack�
9random_rotation_1/rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_3/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_3/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_3StridedSlice-random_rotation_1/rotation_matrix/truediv:z:0@random_rotation_1/rotation_matrix/strided_slice_3/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_3/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_3�
'random_rotation_1/rotation_matrix/Sin_3Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_3�
7random_rotation_1/rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_4/stack�
9random_rotation_1/rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_4/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_4/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_4StridedSlice+random_rotation_1/rotation_matrix/Sin_3:y:0@random_rotation_1/rotation_matrix/strided_slice_4/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_4/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_4�
'random_rotation_1/rotation_matrix/Cos_3Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_3�
7random_rotation_1/rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_5/stack�
9random_rotation_1/rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_5/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_5/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_5StridedSlice+random_rotation_1/rotation_matrix/Cos_3:y:0@random_rotation_1/rotation_matrix/strided_slice_5/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_5/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_5�
7random_rotation_1/rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_6/stack�
9random_rotation_1/rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_6/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_6/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_6StridedSlice/random_rotation_1/rotation_matrix/truediv_1:z:0@random_rotation_1/rotation_matrix/strided_slice_6/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_6/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_6�
-random_rotation_1/rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2/
-random_rotation_1/rotation_matrix/zeros/mul/y�
+random_rotation_1/rotation_matrix/zeros/mulMul8random_rotation_1/rotation_matrix/strided_slice:output:06random_rotation_1/rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2-
+random_rotation_1/rotation_matrix/zeros/mul�
.random_rotation_1/rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�20
.random_rotation_1/rotation_matrix/zeros/Less/y�
,random_rotation_1/rotation_matrix/zeros/LessLess/random_rotation_1/rotation_matrix/zeros/mul:z:07random_rotation_1/rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2.
,random_rotation_1/rotation_matrix/zeros/Less�
0random_rotation_1/rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :22
0random_rotation_1/rotation_matrix/zeros/packed/1�
.random_rotation_1/rotation_matrix/zeros/packedPack8random_rotation_1/rotation_matrix/strided_slice:output:09random_rotation_1/rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:20
.random_rotation_1/rotation_matrix/zeros/packed�
-random_rotation_1/rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2/
-random_rotation_1/rotation_matrix/zeros/Const�
'random_rotation_1/rotation_matrix/zerosFill7random_rotation_1/rotation_matrix/zeros/packed:output:06random_rotation_1/rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/zeros�
-random_rotation_1/rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2/
-random_rotation_1/rotation_matrix/concat/axis�
(random_rotation_1/rotation_matrix/concatConcatV2:random_rotation_1/rotation_matrix/strided_slice_1:output:0)random_rotation_1/rotation_matrix/Neg:y:0:random_rotation_1/rotation_matrix/strided_slice_3:output:0:random_rotation_1/rotation_matrix/strided_slice_4:output:0:random_rotation_1/rotation_matrix/strided_slice_5:output:0:random_rotation_1/rotation_matrix/strided_slice_6:output:00random_rotation_1/rotation_matrix/zeros:output:06random_rotation_1/rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2*
(random_rotation_1/rotation_matrix/concat�
!random_rotation_1/transform/ShapeShape,random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:2#
!random_rotation_1/transform/Shape�
/random_rotation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:21
/random_rotation_1/transform/strided_slice/stack�
1random_rotation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:23
1random_rotation_1/transform/strided_slice/stack_1�
1random_rotation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:23
1random_rotation_1/transform/strided_slice/stack_2�
)random_rotation_1/transform/strided_sliceStridedSlice*random_rotation_1/transform/Shape:output:08random_rotation_1/transform/strided_slice/stack:output:0:random_rotation_1/transform/strided_slice/stack_1:output:0:random_rotation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2+
)random_rotation_1/transform/strided_slice�
&random_rotation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2(
&random_rotation_1/transform/fill_value�
6random_rotation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3,random_flip_1/random_flip_left_right/add:z:01random_rotation_1/rotation_matrix/concat:output:02random_rotation_1/transform/strided_slice:output:0/random_rotation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR28
6random_rotation_1/transform/ImageProjectiveTransformV3�
random_translation_1/ShapeShapeKrandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_translation_1/Shape�
(random_translation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2*
(random_translation_1/strided_slice/stack�
*random_translation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice/stack_1�
*random_translation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice/stack_2�
"random_translation_1/strided_sliceStridedSlice#random_translation_1/Shape:output:01random_translation_1/strided_slice/stack:output:03random_translation_1/strided_slice/stack_1:output:03random_translation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"random_translation_1/strided_slice�
*random_translation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice_1/stack�
,random_translation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_1/stack_1�
,random_translation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_1/stack_2�
$random_translation_1/strided_slice_1StridedSlice#random_translation_1/Shape:output:03random_translation_1/strided_slice_1/stack:output:05random_translation_1/strided_slice_1/stack_1:output:05random_translation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2&
$random_translation_1/strided_slice_1�
random_translation_1/CastCast-random_translation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_translation_1/Cast�
*random_translation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice_2/stack�
,random_translation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_2/stack_1�
,random_translation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_2/stack_2�
$random_translation_1/strided_slice_2StridedSlice#random_translation_1/Shape:output:03random_translation_1/strided_slice_2/stack:output:05random_translation_1/strided_slice_2/stack_1:output:05random_translation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2&
$random_translation_1/strided_slice_2�
random_translation_1/Cast_1Cast-random_translation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_translation_1/Cast_1�
-random_translation_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2/
-random_translation_1/stateful_uniform/shape/1�
+random_translation_1/stateful_uniform/shapePack+random_translation_1/strided_slice:output:06random_translation_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2-
+random_translation_1/stateful_uniform/shape�
)random_translation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2+
)random_translation_1/stateful_uniform/min�
)random_translation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *��>2+
)random_translation_1/stateful_uniform/max�
?random_translation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2A
?random_translation_1/stateful_uniform/StatefulUniform/algorithm�
5random_translation_1/stateful_uniform/StatefulUniformStatefulUniform>random_translation_1_stateful_uniform_statefuluniform_resourceHrandom_translation_1/stateful_uniform/StatefulUniform/algorithm:output:04random_translation_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype027
5random_translation_1/stateful_uniform/StatefulUniform�
)random_translation_1/stateful_uniform/subSub2random_translation_1/stateful_uniform/max:output:02random_translation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2+
)random_translation_1/stateful_uniform/sub�
)random_translation_1/stateful_uniform/mulMul>random_translation_1/stateful_uniform/StatefulUniform:output:0-random_translation_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������2+
)random_translation_1/stateful_uniform/mul�
%random_translation_1/stateful_uniformAdd-random_translation_1/stateful_uniform/mul:z:02random_translation_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2'
%random_translation_1/stateful_uniform�
random_translation_1/mulMul)random_translation_1/stateful_uniform:z:0random_translation_1/Cast:y:0*
T0*'
_output_shapes
:���������2
random_translation_1/mul�
/random_translation_1/stateful_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :21
/random_translation_1/stateful_uniform_1/shape/1�
-random_translation_1/stateful_uniform_1/shapePack+random_translation_1/strided_slice:output:08random_translation_1/stateful_uniform_1/shape/1:output:0*
N*
T0*
_output_shapes
:2/
-random_translation_1/stateful_uniform_1/shape�
+random_translation_1/stateful_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2-
+random_translation_1/stateful_uniform_1/min�
+random_translation_1/stateful_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2-
+random_translation_1/stateful_uniform_1/max�
Arandom_translation_1/stateful_uniform_1/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2C
Arandom_translation_1/stateful_uniform_1/StatefulUniform/algorithm�
7random_translation_1/stateful_uniform_1/StatefulUniformStatefulUniform>random_translation_1_stateful_uniform_statefuluniform_resourceJrandom_translation_1/stateful_uniform_1/StatefulUniform/algorithm:output:06random_translation_1/stateful_uniform_1/shape:output:06^random_translation_1/stateful_uniform/StatefulUniform*'
_output_shapes
:���������*
shape_dtype029
7random_translation_1/stateful_uniform_1/StatefulUniform�
+random_translation_1/stateful_uniform_1/subSub4random_translation_1/stateful_uniform_1/max:output:04random_translation_1/stateful_uniform_1/min:output:0*
T0*
_output_shapes
: 2-
+random_translation_1/stateful_uniform_1/sub�
+random_translation_1/stateful_uniform_1/mulMul@random_translation_1/stateful_uniform_1/StatefulUniform:output:0/random_translation_1/stateful_uniform_1/sub:z:0*
T0*'
_output_shapes
:���������2-
+random_translation_1/stateful_uniform_1/mul�
'random_translation_1/stateful_uniform_1Add/random_translation_1/stateful_uniform_1/mul:z:04random_translation_1/stateful_uniform_1/min:output:0*
T0*'
_output_shapes
:���������2)
'random_translation_1/stateful_uniform_1�
random_translation_1/mul_1Mul+random_translation_1/stateful_uniform_1:z:0random_translation_1/Cast_1:y:0*
T0*'
_output_shapes
:���������2
random_translation_1/mul_1�
 random_translation_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2"
 random_translation_1/concat/axis�
random_translation_1/concatConcatV2random_translation_1/mul_1:z:0random_translation_1/mul:z:0)random_translation_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
random_translation_1/concat�
-random_translation_1/translation_matrix/ShapeShape$random_translation_1/concat:output:0*
T0*
_output_shapes
:2/
-random_translation_1/translation_matrix/Shape�
;random_translation_1/translation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2=
;random_translation_1/translation_matrix/strided_slice/stack�
=random_translation_1/translation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2?
=random_translation_1/translation_matrix/strided_slice/stack_1�
=random_translation_1/translation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2?
=random_translation_1/translation_matrix/strided_slice/stack_2�
5random_translation_1/translation_matrix/strided_sliceStridedSlice6random_translation_1/translation_matrix/Shape:output:0Drandom_translation_1/translation_matrix/strided_slice/stack:output:0Frandom_translation_1/translation_matrix/strided_slice/stack_1:output:0Frandom_translation_1/translation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask27
5random_translation_1/translation_matrix/strided_slice�
2random_translation_1/translation_matrix/ones/mul/yConst*
_output_shapes
: *
dtype0*
value	B :24
2random_translation_1/translation_matrix/ones/mul/y�
0random_translation_1/translation_matrix/ones/mulMul>random_translation_1/translation_matrix/strided_slice:output:0;random_translation_1/translation_matrix/ones/mul/y:output:0*
T0*
_output_shapes
: 22
0random_translation_1/translation_matrix/ones/mul�
3random_translation_1/translation_matrix/ones/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�25
3random_translation_1/translation_matrix/ones/Less/y�
1random_translation_1/translation_matrix/ones/LessLess4random_translation_1/translation_matrix/ones/mul:z:0<random_translation_1/translation_matrix/ones/Less/y:output:0*
T0*
_output_shapes
: 23
1random_translation_1/translation_matrix/ones/Less�
5random_translation_1/translation_matrix/ones/packed/1Const*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/ones/packed/1�
3random_translation_1/translation_matrix/ones/packedPack>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/ones/packed/1:output:0*
N*
T0*
_output_shapes
:25
3random_translation_1/translation_matrix/ones/packed�
2random_translation_1/translation_matrix/ones/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?24
2random_translation_1/translation_matrix/ones/Const�
,random_translation_1/translation_matrix/onesFill<random_translation_1/translation_matrix/ones/packed:output:0;random_translation_1/translation_matrix/ones/Const:output:0*
T0*'
_output_shapes
:���������2.
,random_translation_1/translation_matrix/ones�
3random_translation_1/translation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :25
3random_translation_1/translation_matrix/zeros/mul/y�
1random_translation_1/translation_matrix/zeros/mulMul>random_translation_1/translation_matrix/strided_slice:output:0<random_translation_1/translation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 23
1random_translation_1/translation_matrix/zeros/mul�
4random_translation_1/translation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�26
4random_translation_1/translation_matrix/zeros/Less/y�
2random_translation_1/translation_matrix/zeros/LessLess5random_translation_1/translation_matrix/zeros/mul:z:0=random_translation_1/translation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 24
2random_translation_1/translation_matrix/zeros/Less�
6random_translation_1/translation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :28
6random_translation_1/translation_matrix/zeros/packed/1�
4random_translation_1/translation_matrix/zeros/packedPack>random_translation_1/translation_matrix/strided_slice:output:0?random_translation_1/translation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:26
4random_translation_1/translation_matrix/zeros/packed�
3random_translation_1/translation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    25
3random_translation_1/translation_matrix/zeros/Const�
-random_translation_1/translation_matrix/zerosFill=random_translation_1/translation_matrix/zeros/packed:output:0<random_translation_1/translation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2/
-random_translation_1/translation_matrix/zeros�
=random_translation_1/translation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2?
=random_translation_1/translation_matrix/strided_slice_1/stack�
?random_translation_1/translation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2A
?random_translation_1/translation_matrix/strided_slice_1/stack_1�
?random_translation_1/translation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2A
?random_translation_1/translation_matrix/strided_slice_1/stack_2�
7random_translation_1/translation_matrix/strided_slice_1StridedSlice$random_translation_1/concat:output:0Frandom_translation_1/translation_matrix/strided_slice_1/stack:output:0Hrandom_translation_1/translation_matrix/strided_slice_1/stack_1:output:0Hrandom_translation_1/translation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask29
7random_translation_1/translation_matrix/strided_slice_1�
+random_translation_1/translation_matrix/NegNeg@random_translation_1/translation_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2-
+random_translation_1/translation_matrix/Neg�
5random_translation_1/translation_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/zeros_1/mul/y�
3random_translation_1/translation_matrix/zeros_1/mulMul>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/zeros_1/mul�
6random_translation_1/translation_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�28
6random_translation_1/translation_matrix/zeros_1/Less/y�
4random_translation_1/translation_matrix/zeros_1/LessLess7random_translation_1/translation_matrix/zeros_1/mul:z:0?random_translation_1/translation_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 26
4random_translation_1/translation_matrix/zeros_1/Less�
8random_translation_1/translation_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2:
8random_translation_1/translation_matrix/zeros_1/packed/1�
6random_translation_1/translation_matrix/zeros_1/packedPack>random_translation_1/translation_matrix/strided_slice:output:0Arandom_translation_1/translation_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:28
6random_translation_1/translation_matrix/zeros_1/packed�
5random_translation_1/translation_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    27
5random_translation_1/translation_matrix/zeros_1/Const�
/random_translation_1/translation_matrix/zeros_1Fill?random_translation_1/translation_matrix/zeros_1/packed:output:0>random_translation_1/translation_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������21
/random_translation_1/translation_matrix/zeros_1�
4random_translation_1/translation_matrix/ones_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :26
4random_translation_1/translation_matrix/ones_1/mul/y�
2random_translation_1/translation_matrix/ones_1/mulMul>random_translation_1/translation_matrix/strided_slice:output:0=random_translation_1/translation_matrix/ones_1/mul/y:output:0*
T0*
_output_shapes
: 24
2random_translation_1/translation_matrix/ones_1/mul�
5random_translation_1/translation_matrix/ones_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�27
5random_translation_1/translation_matrix/ones_1/Less/y�
3random_translation_1/translation_matrix/ones_1/LessLess6random_translation_1/translation_matrix/ones_1/mul:z:0>random_translation_1/translation_matrix/ones_1/Less/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/ones_1/Less�
7random_translation_1/translation_matrix/ones_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :29
7random_translation_1/translation_matrix/ones_1/packed/1�
5random_translation_1/translation_matrix/ones_1/packedPack>random_translation_1/translation_matrix/strided_slice:output:0@random_translation_1/translation_matrix/ones_1/packed/1:output:0*
N*
T0*
_output_shapes
:27
5random_translation_1/translation_matrix/ones_1/packed�
4random_translation_1/translation_matrix/ones_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?26
4random_translation_1/translation_matrix/ones_1/Const�
.random_translation_1/translation_matrix/ones_1Fill>random_translation_1/translation_matrix/ones_1/packed:output:0=random_translation_1/translation_matrix/ones_1/Const:output:0*
T0*'
_output_shapes
:���������20
.random_translation_1/translation_matrix/ones_1�
=random_translation_1/translation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2?
=random_translation_1/translation_matrix/strided_slice_2/stack�
?random_translation_1/translation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2A
?random_translation_1/translation_matrix/strided_slice_2/stack_1�
?random_translation_1/translation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2A
?random_translation_1/translation_matrix/strided_slice_2/stack_2�
7random_translation_1/translation_matrix/strided_slice_2StridedSlice$random_translation_1/concat:output:0Frandom_translation_1/translation_matrix/strided_slice_2/stack:output:0Hrandom_translation_1/translation_matrix/strided_slice_2/stack_1:output:0Hrandom_translation_1/translation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask29
7random_translation_1/translation_matrix/strided_slice_2�
-random_translation_1/translation_matrix/Neg_1Neg@random_translation_1/translation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2/
-random_translation_1/translation_matrix/Neg_1�
5random_translation_1/translation_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/zeros_2/mul/y�
3random_translation_1/translation_matrix/zeros_2/mulMul>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/zeros_2/mul�
6random_translation_1/translation_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�28
6random_translation_1/translation_matrix/zeros_2/Less/y�
4random_translation_1/translation_matrix/zeros_2/LessLess7random_translation_1/translation_matrix/zeros_2/mul:z:0?random_translation_1/translation_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 26
4random_translation_1/translation_matrix/zeros_2/Less�
8random_translation_1/translation_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2:
8random_translation_1/translation_matrix/zeros_2/packed/1�
6random_translation_1/translation_matrix/zeros_2/packedPack>random_translation_1/translation_matrix/strided_slice:output:0Arandom_translation_1/translation_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:28
6random_translation_1/translation_matrix/zeros_2/packed�
5random_translation_1/translation_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    27
5random_translation_1/translation_matrix/zeros_2/Const�
/random_translation_1/translation_matrix/zeros_2Fill?random_translation_1/translation_matrix/zeros_2/packed:output:0>random_translation_1/translation_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������21
/random_translation_1/translation_matrix/zeros_2�
3random_translation_1/translation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :25
3random_translation_1/translation_matrix/concat/axis�
.random_translation_1/translation_matrix/concatConcatV25random_translation_1/translation_matrix/ones:output:06random_translation_1/translation_matrix/zeros:output:0/random_translation_1/translation_matrix/Neg:y:08random_translation_1/translation_matrix/zeros_1:output:07random_translation_1/translation_matrix/ones_1:output:01random_translation_1/translation_matrix/Neg_1:y:08random_translation_1/translation_matrix/zeros_2:output:0<random_translation_1/translation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������20
.random_translation_1/translation_matrix/concat�
$random_translation_1/transform/ShapeShapeKrandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2&
$random_translation_1/transform/Shape�
2random_translation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:24
2random_translation_1/transform/strided_slice/stack�
4random_translation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:26
4random_translation_1/transform/strided_slice/stack_1�
4random_translation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:26
4random_translation_1/transform/strided_slice/stack_2�
,random_translation_1/transform/strided_sliceStridedSlice-random_translation_1/transform/Shape:output:0;random_translation_1/transform/strided_slice/stack:output:0=random_translation_1/transform/strided_slice/stack_1:output:0=random_translation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2.
,random_translation_1/transform/strided_slice�
)random_translation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2+
)random_translation_1/transform/fill_value�
9random_translation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Krandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:07random_translation_1/translation_matrix/concat:output:05random_translation_1/transform/strided_slice:output:02random_translation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	NEAREST*
interpolation
BILINEAR2;
9random_translation_1/transform/ImageProjectiveTransformV3�
random_zoom_1/ShapeShapeNrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_zoom_1/Shape�
!random_zoom_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2#
!random_zoom_1/strided_slice/stack�
#random_zoom_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice/stack_1�
#random_zoom_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice/stack_2�
random_zoom_1/strided_sliceStridedSlicerandom_zoom_1/Shape:output:0*random_zoom_1/strided_slice/stack:output:0,random_zoom_1/strided_slice/stack_1:output:0,random_zoom_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice�
#random_zoom_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice_1/stack�
%random_zoom_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_1/stack_1�
%random_zoom_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_1/stack_2�
random_zoom_1/strided_slice_1StridedSlicerandom_zoom_1/Shape:output:0,random_zoom_1/strided_slice_1/stack:output:0.random_zoom_1/strided_slice_1/stack_1:output:0.random_zoom_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice_1�
random_zoom_1/CastCast&random_zoom_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_zoom_1/Cast�
#random_zoom_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice_2/stack�
%random_zoom_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_2/stack_1�
%random_zoom_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_2/stack_2�
random_zoom_1/strided_slice_2StridedSlicerandom_zoom_1/Shape:output:0,random_zoom_1/strided_slice_2/stack:output:0.random_zoom_1/strided_slice_2/stack_1:output:0.random_zoom_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice_2�
random_zoom_1/Cast_1Cast&random_zoom_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_zoom_1/Cast_1�
&random_zoom_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2(
&random_zoom_1/stateful_uniform/shape/1�
$random_zoom_1/stateful_uniform/shapePack$random_zoom_1/strided_slice:output:0/random_zoom_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2&
$random_zoom_1/stateful_uniform/shape�
"random_zoom_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *��L?2$
"random_zoom_1/stateful_uniform/min�
"random_zoom_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���?2$
"random_zoom_1/stateful_uniform/max�
8random_zoom_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2:
8random_zoom_1/stateful_uniform/StatefulUniform/algorithm�
.random_zoom_1/stateful_uniform/StatefulUniformStatefulUniform7random_zoom_1_stateful_uniform_statefuluniform_resourceArandom_zoom_1/stateful_uniform/StatefulUniform/algorithm:output:0-random_zoom_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype020
.random_zoom_1/stateful_uniform/StatefulUniform�
"random_zoom_1/stateful_uniform/subSub+random_zoom_1/stateful_uniform/max:output:0+random_zoom_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2$
"random_zoom_1/stateful_uniform/sub�
"random_zoom_1/stateful_uniform/mulMul7random_zoom_1/stateful_uniform/StatefulUniform:output:0&random_zoom_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������2$
"random_zoom_1/stateful_uniform/mul�
random_zoom_1/stateful_uniformAdd&random_zoom_1/stateful_uniform/mul:z:0+random_zoom_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2 
random_zoom_1/stateful_uniformx
random_zoom_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
random_zoom_1/concat/axis�
random_zoom_1/concatConcatV2"random_zoom_1/stateful_uniform:z:0"random_zoom_1/stateful_uniform:z:0"random_zoom_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
random_zoom_1/concat�
random_zoom_1/zoom_matrix/ShapeShaperandom_zoom_1/concat:output:0*
T0*
_output_shapes
:2!
random_zoom_1/zoom_matrix/Shape�
-random_zoom_1/zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2/
-random_zoom_1/zoom_matrix/strided_slice/stack�
/random_zoom_1/zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:21
/random_zoom_1/zoom_matrix/strided_slice/stack_1�
/random_zoom_1/zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:21
/random_zoom_1/zoom_matrix/strided_slice/stack_2�
'random_zoom_1/zoom_matrix/strided_sliceStridedSlice(random_zoom_1/zoom_matrix/Shape:output:06random_zoom_1/zoom_matrix/strided_slice/stack:output:08random_zoom_1/zoom_matrix/strided_slice/stack_1:output:08random_zoom_1/zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2)
'random_zoom_1/zoom_matrix/strided_slice�
random_zoom_1/zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2!
random_zoom_1/zoom_matrix/sub/y�
random_zoom_1/zoom_matrix/subSubrandom_zoom_1/Cast_1:y:0(random_zoom_1/zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 2
random_zoom_1/zoom_matrix/sub�
#random_zoom_1/zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2%
#random_zoom_1/zoom_matrix/truediv/y�
!random_zoom_1/zoom_matrix/truedivRealDiv!random_zoom_1/zoom_matrix/sub:z:0,random_zoom_1/zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 2#
!random_zoom_1/zoom_matrix/truediv�
/random_zoom_1/zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            21
/random_zoom_1/zoom_matrix/strided_slice_1/stack�
1random_zoom_1/zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_1/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_1/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_1StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_1/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_1/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_1�
!random_zoom_1/zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_1/x�
random_zoom_1/zoom_matrix/sub_1Sub*random_zoom_1/zoom_matrix/sub_1/x:output:02random_zoom_1/zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/sub_1�
random_zoom_1/zoom_matrix/mulMul%random_zoom_1/zoom_matrix/truediv:z:0#random_zoom_1/zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:���������2
random_zoom_1/zoom_matrix/mul�
!random_zoom_1/zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_2/y�
random_zoom_1/zoom_matrix/sub_2Subrandom_zoom_1/Cast:y:0*random_zoom_1/zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2!
random_zoom_1/zoom_matrix/sub_2�
%random_zoom_1/zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2'
%random_zoom_1/zoom_matrix/truediv_1/y�
#random_zoom_1/zoom_matrix/truediv_1RealDiv#random_zoom_1/zoom_matrix/sub_2:z:0.random_zoom_1/zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 2%
#random_zoom_1/zoom_matrix/truediv_1�
/random_zoom_1/zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom_1/zoom_matrix/strided_slice_2/stack�
1random_zoom_1/zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_2/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_2/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_2StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_2/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_2/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_2�
!random_zoom_1/zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_3/x�
random_zoom_1/zoom_matrix/sub_3Sub*random_zoom_1/zoom_matrix/sub_3/x:output:02random_zoom_1/zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/sub_3�
random_zoom_1/zoom_matrix/mul_1Mul'random_zoom_1/zoom_matrix/truediv_1:z:0#random_zoom_1/zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/mul_1�
/random_zoom_1/zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            21
/random_zoom_1/zoom_matrix/strided_slice_3/stack�
1random_zoom_1/zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_3/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_3/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_3StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_3/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_3/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_3�
%random_zoom_1/zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2'
%random_zoom_1/zoom_matrix/zeros/mul/y�
#random_zoom_1/zoom_matrix/zeros/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:0.random_zoom_1/zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2%
#random_zoom_1/zoom_matrix/zeros/mul�
&random_zoom_1/zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2(
&random_zoom_1/zoom_matrix/zeros/Less/y�
$random_zoom_1/zoom_matrix/zeros/LessLess'random_zoom_1/zoom_matrix/zeros/mul:z:0/random_zoom_1/zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2&
$random_zoom_1/zoom_matrix/zeros/Less�
(random_zoom_1/zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2*
(random_zoom_1/zoom_matrix/zeros/packed/1�
&random_zoom_1/zoom_matrix/zeros/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:01random_zoom_1/zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2(
&random_zoom_1/zoom_matrix/zeros/packed�
%random_zoom_1/zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%random_zoom_1/zoom_matrix/zeros/Const�
random_zoom_1/zoom_matrix/zerosFill/random_zoom_1/zoom_matrix/zeros/packed:output:0.random_zoom_1/zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/zeros�
'random_zoom_1/zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2)
'random_zoom_1/zoom_matrix/zeros_1/mul/y�
%random_zoom_1/zoom_matrix/zeros_1/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:00random_zoom_1/zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2'
%random_zoom_1/zoom_matrix/zeros_1/mul�
(random_zoom_1/zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2*
(random_zoom_1/zoom_matrix/zeros_1/Less/y�
&random_zoom_1/zoom_matrix/zeros_1/LessLess)random_zoom_1/zoom_matrix/zeros_1/mul:z:01random_zoom_1/zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2(
&random_zoom_1/zoom_matrix/zeros_1/Less�
*random_zoom_1/zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2,
*random_zoom_1/zoom_matrix/zeros_1/packed/1�
(random_zoom_1/zoom_matrix/zeros_1/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:03random_zoom_1/zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2*
(random_zoom_1/zoom_matrix/zeros_1/packed�
'random_zoom_1/zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2)
'random_zoom_1/zoom_matrix/zeros_1/Const�
!random_zoom_1/zoom_matrix/zeros_1Fill1random_zoom_1/zoom_matrix/zeros_1/packed:output:00random_zoom_1/zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������2#
!random_zoom_1/zoom_matrix/zeros_1�
/random_zoom_1/zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom_1/zoom_matrix/strided_slice_4/stack�
1random_zoom_1/zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_4/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_4/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_4StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_4/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_4/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_4�
'random_zoom_1/zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2)
'random_zoom_1/zoom_matrix/zeros_2/mul/y�
%random_zoom_1/zoom_matrix/zeros_2/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:00random_zoom_1/zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2'
%random_zoom_1/zoom_matrix/zeros_2/mul�
(random_zoom_1/zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2*
(random_zoom_1/zoom_matrix/zeros_2/Less/y�
&random_zoom_1/zoom_matrix/zeros_2/LessLess)random_zoom_1/zoom_matrix/zeros_2/mul:z:01random_zoom_1/zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2(
&random_zoom_1/zoom_matrix/zeros_2/Less�
*random_zoom_1/zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2,
*random_zoom_1/zoom_matrix/zeros_2/packed/1�
(random_zoom_1/zoom_matrix/zeros_2/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:03random_zoom_1/zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2*
(random_zoom_1/zoom_matrix/zeros_2/packed�
'random_zoom_1/zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2)
'random_zoom_1/zoom_matrix/zeros_2/Const�
!random_zoom_1/zoom_matrix/zeros_2Fill1random_zoom_1/zoom_matrix/zeros_2/packed:output:00random_zoom_1/zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������2#
!random_zoom_1/zoom_matrix/zeros_2�
%random_zoom_1/zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2'
%random_zoom_1/zoom_matrix/concat/axis�
 random_zoom_1/zoom_matrix/concatConcatV22random_zoom_1/zoom_matrix/strided_slice_3:output:0(random_zoom_1/zoom_matrix/zeros:output:0!random_zoom_1/zoom_matrix/mul:z:0*random_zoom_1/zoom_matrix/zeros_1:output:02random_zoom_1/zoom_matrix/strided_slice_4:output:0#random_zoom_1/zoom_matrix/mul_1:z:0*random_zoom_1/zoom_matrix/zeros_2:output:0.random_zoom_1/zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2"
 random_zoom_1/zoom_matrix/concat�
random_zoom_1/transform/ShapeShapeNrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_zoom_1/transform/Shape�
+random_zoom_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2-
+random_zoom_1/transform/strided_slice/stack�
-random_zoom_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2/
-random_zoom_1/transform/strided_slice/stack_1�
-random_zoom_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2/
-random_zoom_1/transform/strided_slice/stack_2�
%random_zoom_1/transform/strided_sliceStridedSlice&random_zoom_1/transform/Shape:output:04random_zoom_1/transform/strided_slice/stack:output:06random_zoom_1/transform/strided_slice/stack_1:output:06random_zoom_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2'
%random_zoom_1/transform/strided_slice�
"random_zoom_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2$
"random_zoom_1/transform/fill_value�
2random_zoom_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Nrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0)random_zoom_1/zoom_matrix/concat:output:0.random_zoom_1/transform/strided_slice:output:0+random_zoom_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR24
2random_zoom_1/transform/ImageProjectiveTransformV3�
IdentityIdentityGrandom_zoom_1/transform/ImageProjectiveTransformV3:transformed_images:03^random_rotation_1/stateful_uniform/StatefulUniform6^random_translation_1/stateful_uniform/StatefulUniform8^random_translation_1/stateful_uniform_1/StatefulUniform/^random_zoom_1/stateful_uniform/StatefulUniform*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*:
_input_shapes)
':���������&:::2h
2random_rotation_1/stateful_uniform/StatefulUniform2random_rotation_1/stateful_uniform/StatefulUniform2n
5random_translation_1/stateful_uniform/StatefulUniform5random_translation_1/stateful_uniform/StatefulUniform2r
7random_translation_1/stateful_uniform_1/StatefulUniform7random_translation_1/stateful_uniform_1/StatefulUniform2`
.random_zoom_1/stateful_uniform/StatefulUniform.random_zoom_1/stateful_uniform/StatefulUniform:d `
/
_output_shapes
:���������&
-
_user_specified_namerandom_flip_1_input
�
f
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_31317

inputs
identity�
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4������������������������������������*
ksize
*
paddingVALID*
strides
2	
MaxPool�
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4������������������������������������2

Identity"
identityIdentity:output:0*I
_input_shapes8
6:4������������������������������������:r n
J
_output_shapes8
6:4������������������������������������
 
_user_specified_nameinputs
��
�
G__inference_sequential_3_layer_call_and_return_conditional_losses_31293

inputs?
;random_rotation_1_stateful_uniform_statefuluniform_resourceB
>random_translation_1_stateful_uniform_statefuluniform_resource;
7random_zoom_1_stateful_uniform_statefuluniform_resource
identity��2random_rotation_1/stateful_uniform/StatefulUniform�5random_translation_1/stateful_uniform/StatefulUniform�7random_translation_1/stateful_uniform_1/StatefulUniform�.random_zoom_1/stateful_uniform/StatefulUniform�
7random_flip_1/random_flip_left_right/control_dependencyIdentityinputs*
T0*
_class
loc:@inputs*/
_output_shapes
:���������&29
7random_flip_1/random_flip_left_right/control_dependency�
*random_flip_1/random_flip_left_right/ShapeShape@random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:2,
*random_flip_1/random_flip_left_right/Shape�
8random_flip_1/random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2:
8random_flip_1/random_flip_left_right/strided_slice/stack�
:random_flip_1/random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2<
:random_flip_1/random_flip_left_right/strided_slice/stack_1�
:random_flip_1/random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2<
:random_flip_1/random_flip_left_right/strided_slice/stack_2�
2random_flip_1/random_flip_left_right/strided_sliceStridedSlice3random_flip_1/random_flip_left_right/Shape:output:0Arandom_flip_1/random_flip_left_right/strided_slice/stack:output:0Crandom_flip_1/random_flip_left_right/strided_slice/stack_1:output:0Crandom_flip_1/random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask24
2random_flip_1/random_flip_left_right/strided_slice�
9random_flip_1/random_flip_left_right/random_uniform/shapePack;random_flip_1/random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2;
9random_flip_1/random_flip_left_right/random_uniform/shape�
7random_flip_1/random_flip_left_right/random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    29
7random_flip_1/random_flip_left_right/random_uniform/min�
7random_flip_1/random_flip_left_right/random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  �?29
7random_flip_1/random_flip_left_right/random_uniform/max�
Arandom_flip_1/random_flip_left_right/random_uniform/RandomUniformRandomUniformBrandom_flip_1/random_flip_left_right/random_uniform/shape:output:0*
T0*#
_output_shapes
:���������*
dtype02C
Arandom_flip_1/random_flip_left_right/random_uniform/RandomUniform�
7random_flip_1/random_flip_left_right/random_uniform/MulMulJrandom_flip_1/random_flip_left_right/random_uniform/RandomUniform:output:0@random_flip_1/random_flip_left_right/random_uniform/max:output:0*
T0*#
_output_shapes
:���������29
7random_flip_1/random_flip_left_right/random_uniform/Mul�
4random_flip_1/random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/1�
4random_flip_1/random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/2�
4random_flip_1/random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :26
4random_flip_1/random_flip_left_right/Reshape/shape/3�
2random_flip_1/random_flip_left_right/Reshape/shapePack;random_flip_1/random_flip_left_right/strided_slice:output:0=random_flip_1/random_flip_left_right/Reshape/shape/1:output:0=random_flip_1/random_flip_left_right/Reshape/shape/2:output:0=random_flip_1/random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:24
2random_flip_1/random_flip_left_right/Reshape/shape�
,random_flip_1/random_flip_left_right/ReshapeReshape;random_flip_1/random_flip_left_right/random_uniform/Mul:z:0;random_flip_1/random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:���������2.
,random_flip_1/random_flip_left_right/Reshape�
*random_flip_1/random_flip_left_right/RoundRound5random_flip_1/random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:���������2,
*random_flip_1/random_flip_left_right/Round�
3random_flip_1/random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:25
3random_flip_1/random_flip_left_right/ReverseV2/axis�
.random_flip_1/random_flip_left_right/ReverseV2	ReverseV2@random_flip_1/random_flip_left_right/control_dependency:output:0<random_flip_1/random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:���������&20
.random_flip_1/random_flip_left_right/ReverseV2�
(random_flip_1/random_flip_left_right/mulMul.random_flip_1/random_flip_left_right/Round:y:07random_flip_1/random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:���������&2*
(random_flip_1/random_flip_left_right/mul�
*random_flip_1/random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2,
*random_flip_1/random_flip_left_right/sub/x�
(random_flip_1/random_flip_left_right/subSub3random_flip_1/random_flip_left_right/sub/x:output:0.random_flip_1/random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:���������2*
(random_flip_1/random_flip_left_right/sub�
*random_flip_1/random_flip_left_right/mul_1Mul,random_flip_1/random_flip_left_right/sub:z:0@random_flip_1/random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:���������&2,
*random_flip_1/random_flip_left_right/mul_1�
(random_flip_1/random_flip_left_right/addAddV2,random_flip_1/random_flip_left_right/mul:z:0.random_flip_1/random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:���������&2*
(random_flip_1/random_flip_left_right/add�
random_rotation_1/ShapeShape,random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:2
random_rotation_1/Shape�
%random_rotation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2'
%random_rotation_1/strided_slice/stack�
'random_rotation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice/stack_1�
'random_rotation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice/stack_2�
random_rotation_1/strided_sliceStridedSlice random_rotation_1/Shape:output:0.random_rotation_1/strided_slice/stack:output:00random_rotation_1/strided_slice/stack_1:output:00random_rotation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2!
random_rotation_1/strided_slice�
'random_rotation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice_1/stack�
)random_rotation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_1/stack_1�
)random_rotation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_1/stack_2�
!random_rotation_1/strided_slice_1StridedSlice random_rotation_1/Shape:output:00random_rotation_1/strided_slice_1/stack:output:02random_rotation_1/strided_slice_1/stack_1:output:02random_rotation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2#
!random_rotation_1/strided_slice_1�
random_rotation_1/CastCast*random_rotation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_rotation_1/Cast�
'random_rotation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation_1/strided_slice_2/stack�
)random_rotation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_2/stack_1�
)random_rotation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2+
)random_rotation_1/strided_slice_2/stack_2�
!random_rotation_1/strided_slice_2StridedSlice random_rotation_1/Shape:output:00random_rotation_1/strided_slice_2/stack:output:02random_rotation_1/strided_slice_2/stack_1:output:02random_rotation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2#
!random_rotation_1/strided_slice_2�
random_rotation_1/Cast_1Cast*random_rotation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_rotation_1/Cast_1�
(random_rotation_1/stateful_uniform/shapePack(random_rotation_1/strided_slice:output:0*
N*
T0*
_output_shapes
:2*
(random_rotation_1/stateful_uniform/shape�
&random_rotation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *   �2(
&random_rotation_1/stateful_uniform/min�
&random_rotation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2(
&random_rotation_1/stateful_uniform/max�
<random_rotation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2>
<random_rotation_1/stateful_uniform/StatefulUniform/algorithm�
2random_rotation_1/stateful_uniform/StatefulUniformStatefulUniform;random_rotation_1_stateful_uniform_statefuluniform_resourceErandom_rotation_1/stateful_uniform/StatefulUniform/algorithm:output:01random_rotation_1/stateful_uniform/shape:output:0*#
_output_shapes
:���������*
shape_dtype024
2random_rotation_1/stateful_uniform/StatefulUniform�
&random_rotation_1/stateful_uniform/subSub/random_rotation_1/stateful_uniform/max:output:0/random_rotation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2(
&random_rotation_1/stateful_uniform/sub�
&random_rotation_1/stateful_uniform/mulMul;random_rotation_1/stateful_uniform/StatefulUniform:output:0*random_rotation_1/stateful_uniform/sub:z:0*
T0*#
_output_shapes
:���������2(
&random_rotation_1/stateful_uniform/mul�
"random_rotation_1/stateful_uniformAdd*random_rotation_1/stateful_uniform/mul:z:0/random_rotation_1/stateful_uniform/min:output:0*
T0*#
_output_shapes
:���������2$
"random_rotation_1/stateful_uniform�
'random_rotation_1/rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2)
'random_rotation_1/rotation_matrix/sub/y�
%random_rotation_1/rotation_matrix/subSubrandom_rotation_1/Cast_1:y:00random_rotation_1/rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 2'
%random_rotation_1/rotation_matrix/sub�
%random_rotation_1/rotation_matrix/CosCos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Cos�
)random_rotation_1/rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_1/y�
'random_rotation_1/rotation_matrix/sub_1Subrandom_rotation_1/Cast_1:y:02random_rotation_1/rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_1�
%random_rotation_1/rotation_matrix/mulMul)random_rotation_1/rotation_matrix/Cos:y:0+random_rotation_1/rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/mul�
%random_rotation_1/rotation_matrix/SinSin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Sin�
)random_rotation_1/rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_2/y�
'random_rotation_1/rotation_matrix/sub_2Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_2�
'random_rotation_1/rotation_matrix/mul_1Mul)random_rotation_1/rotation_matrix/Sin:y:0+random_rotation_1/rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_1�
'random_rotation_1/rotation_matrix/sub_3Sub)random_rotation_1/rotation_matrix/mul:z:0+random_rotation_1/rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_3�
'random_rotation_1/rotation_matrix/sub_4Sub)random_rotation_1/rotation_matrix/sub:z:0+random_rotation_1/rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_4�
+random_rotation_1/rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2-
+random_rotation_1/rotation_matrix/truediv/y�
)random_rotation_1/rotation_matrix/truedivRealDiv+random_rotation_1/rotation_matrix/sub_4:z:04random_rotation_1/rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:���������2+
)random_rotation_1/rotation_matrix/truediv�
)random_rotation_1/rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_5/y�
'random_rotation_1/rotation_matrix/sub_5Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_5�
'random_rotation_1/rotation_matrix/Sin_1Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_1�
)random_rotation_1/rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_6/y�
'random_rotation_1/rotation_matrix/sub_6Subrandom_rotation_1/Cast_1:y:02random_rotation_1/rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_6�
'random_rotation_1/rotation_matrix/mul_2Mul+random_rotation_1/rotation_matrix/Sin_1:y:0+random_rotation_1/rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_2�
'random_rotation_1/rotation_matrix/Cos_1Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_1�
)random_rotation_1/rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2+
)random_rotation_1/rotation_matrix/sub_7/y�
'random_rotation_1/rotation_matrix/sub_7Subrandom_rotation_1/Cast:y:02random_rotation_1/rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 2)
'random_rotation_1/rotation_matrix/sub_7�
'random_rotation_1/rotation_matrix/mul_3Mul+random_rotation_1/rotation_matrix/Cos_1:y:0+random_rotation_1/rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/mul_3�
%random_rotation_1/rotation_matrix/addAddV2+random_rotation_1/rotation_matrix/mul_2:z:0+random_rotation_1/rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/add�
'random_rotation_1/rotation_matrix/sub_8Sub+random_rotation_1/rotation_matrix/sub_5:z:0)random_rotation_1/rotation_matrix/add:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/sub_8�
-random_rotation_1/rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2/
-random_rotation_1/rotation_matrix/truediv_1/y�
+random_rotation_1/rotation_matrix/truediv_1RealDiv+random_rotation_1/rotation_matrix/sub_8:z:06random_rotation_1/rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:���������2-
+random_rotation_1/rotation_matrix/truediv_1�
'random_rotation_1/rotation_matrix/ShapeShape&random_rotation_1/stateful_uniform:z:0*
T0*
_output_shapes
:2)
'random_rotation_1/rotation_matrix/Shape�
5random_rotation_1/rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 27
5random_rotation_1/rotation_matrix/strided_slice/stack�
7random_rotation_1/rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:29
7random_rotation_1/rotation_matrix/strided_slice/stack_1�
7random_rotation_1/rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:29
7random_rotation_1/rotation_matrix/strided_slice/stack_2�
/random_rotation_1/rotation_matrix/strided_sliceStridedSlice0random_rotation_1/rotation_matrix/Shape:output:0>random_rotation_1/rotation_matrix/strided_slice/stack:output:0@random_rotation_1/rotation_matrix/strided_slice/stack_1:output:0@random_rotation_1/rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask21
/random_rotation_1/rotation_matrix/strided_slice�
'random_rotation_1/rotation_matrix/Cos_2Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_2�
7random_rotation_1/rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_1/stack�
9random_rotation_1/rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_1/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_1/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_1StridedSlice+random_rotation_1/rotation_matrix/Cos_2:y:0@random_rotation_1/rotation_matrix/strided_slice_1/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_1/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_1�
'random_rotation_1/rotation_matrix/Sin_2Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_2�
7random_rotation_1/rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_2/stack�
9random_rotation_1/rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_2/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_2/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_2StridedSlice+random_rotation_1/rotation_matrix/Sin_2:y:0@random_rotation_1/rotation_matrix/strided_slice_2/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_2/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_2�
%random_rotation_1/rotation_matrix/NegNeg:random_rotation_1/rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2'
%random_rotation_1/rotation_matrix/Neg�
7random_rotation_1/rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_3/stack�
9random_rotation_1/rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_3/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_3/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_3StridedSlice-random_rotation_1/rotation_matrix/truediv:z:0@random_rotation_1/rotation_matrix/strided_slice_3/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_3/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_3�
'random_rotation_1/rotation_matrix/Sin_3Sin&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Sin_3�
7random_rotation_1/rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_4/stack�
9random_rotation_1/rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_4/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_4/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_4StridedSlice+random_rotation_1/rotation_matrix/Sin_3:y:0@random_rotation_1/rotation_matrix/strided_slice_4/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_4/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_4�
'random_rotation_1/rotation_matrix/Cos_3Cos&random_rotation_1/stateful_uniform:z:0*
T0*#
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/Cos_3�
7random_rotation_1/rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_5/stack�
9random_rotation_1/rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_5/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_5/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_5StridedSlice+random_rotation_1/rotation_matrix/Cos_3:y:0@random_rotation_1/rotation_matrix/strided_slice_5/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_5/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_5�
7random_rotation_1/rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation_1/rotation_matrix/strided_slice_6/stack�
9random_rotation_1/rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2;
9random_rotation_1/rotation_matrix/strided_slice_6/stack_1�
9random_rotation_1/rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2;
9random_rotation_1/rotation_matrix/strided_slice_6/stack_2�
1random_rotation_1/rotation_matrix/strided_slice_6StridedSlice/random_rotation_1/rotation_matrix/truediv_1:z:0@random_rotation_1/rotation_matrix/strided_slice_6/stack:output:0Brandom_rotation_1/rotation_matrix/strided_slice_6/stack_1:output:0Brandom_rotation_1/rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask23
1random_rotation_1/rotation_matrix/strided_slice_6�
-random_rotation_1/rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2/
-random_rotation_1/rotation_matrix/zeros/mul/y�
+random_rotation_1/rotation_matrix/zeros/mulMul8random_rotation_1/rotation_matrix/strided_slice:output:06random_rotation_1/rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2-
+random_rotation_1/rotation_matrix/zeros/mul�
.random_rotation_1/rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�20
.random_rotation_1/rotation_matrix/zeros/Less/y�
,random_rotation_1/rotation_matrix/zeros/LessLess/random_rotation_1/rotation_matrix/zeros/mul:z:07random_rotation_1/rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2.
,random_rotation_1/rotation_matrix/zeros/Less�
0random_rotation_1/rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :22
0random_rotation_1/rotation_matrix/zeros/packed/1�
.random_rotation_1/rotation_matrix/zeros/packedPack8random_rotation_1/rotation_matrix/strided_slice:output:09random_rotation_1/rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:20
.random_rotation_1/rotation_matrix/zeros/packed�
-random_rotation_1/rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2/
-random_rotation_1/rotation_matrix/zeros/Const�
'random_rotation_1/rotation_matrix/zerosFill7random_rotation_1/rotation_matrix/zeros/packed:output:06random_rotation_1/rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2)
'random_rotation_1/rotation_matrix/zeros�
-random_rotation_1/rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2/
-random_rotation_1/rotation_matrix/concat/axis�
(random_rotation_1/rotation_matrix/concatConcatV2:random_rotation_1/rotation_matrix/strided_slice_1:output:0)random_rotation_1/rotation_matrix/Neg:y:0:random_rotation_1/rotation_matrix/strided_slice_3:output:0:random_rotation_1/rotation_matrix/strided_slice_4:output:0:random_rotation_1/rotation_matrix/strided_slice_5:output:0:random_rotation_1/rotation_matrix/strided_slice_6:output:00random_rotation_1/rotation_matrix/zeros:output:06random_rotation_1/rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2*
(random_rotation_1/rotation_matrix/concat�
!random_rotation_1/transform/ShapeShape,random_flip_1/random_flip_left_right/add:z:0*
T0*
_output_shapes
:2#
!random_rotation_1/transform/Shape�
/random_rotation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:21
/random_rotation_1/transform/strided_slice/stack�
1random_rotation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:23
1random_rotation_1/transform/strided_slice/stack_1�
1random_rotation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:23
1random_rotation_1/transform/strided_slice/stack_2�
)random_rotation_1/transform/strided_sliceStridedSlice*random_rotation_1/transform/Shape:output:08random_rotation_1/transform/strided_slice/stack:output:0:random_rotation_1/transform/strided_slice/stack_1:output:0:random_rotation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2+
)random_rotation_1/transform/strided_slice�
&random_rotation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2(
&random_rotation_1/transform/fill_value�
6random_rotation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3,random_flip_1/random_flip_left_right/add:z:01random_rotation_1/rotation_matrix/concat:output:02random_rotation_1/transform/strided_slice:output:0/random_rotation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR28
6random_rotation_1/transform/ImageProjectiveTransformV3�
random_translation_1/ShapeShapeKrandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_translation_1/Shape�
(random_translation_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2*
(random_translation_1/strided_slice/stack�
*random_translation_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice/stack_1�
*random_translation_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice/stack_2�
"random_translation_1/strided_sliceStridedSlice#random_translation_1/Shape:output:01random_translation_1/strided_slice/stack:output:03random_translation_1/strided_slice/stack_1:output:03random_translation_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"random_translation_1/strided_slice�
*random_translation_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice_1/stack�
,random_translation_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_1/stack_1�
,random_translation_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_1/stack_2�
$random_translation_1/strided_slice_1StridedSlice#random_translation_1/Shape:output:03random_translation_1/strided_slice_1/stack:output:05random_translation_1/strided_slice_1/stack_1:output:05random_translation_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2&
$random_translation_1/strided_slice_1�
random_translation_1/CastCast-random_translation_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_translation_1/Cast�
*random_translation_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2,
*random_translation_1/strided_slice_2/stack�
,random_translation_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_2/stack_1�
,random_translation_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,random_translation_1/strided_slice_2/stack_2�
$random_translation_1/strided_slice_2StridedSlice#random_translation_1/Shape:output:03random_translation_1/strided_slice_2/stack:output:05random_translation_1/strided_slice_2/stack_1:output:05random_translation_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2&
$random_translation_1/strided_slice_2�
random_translation_1/Cast_1Cast-random_translation_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_translation_1/Cast_1�
-random_translation_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2/
-random_translation_1/stateful_uniform/shape/1�
+random_translation_1/stateful_uniform/shapePack+random_translation_1/strided_slice:output:06random_translation_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2-
+random_translation_1/stateful_uniform/shape�
)random_translation_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2+
)random_translation_1/stateful_uniform/min�
)random_translation_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *��>2+
)random_translation_1/stateful_uniform/max�
?random_translation_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2A
?random_translation_1/stateful_uniform/StatefulUniform/algorithm�
5random_translation_1/stateful_uniform/StatefulUniformStatefulUniform>random_translation_1_stateful_uniform_statefuluniform_resourceHrandom_translation_1/stateful_uniform/StatefulUniform/algorithm:output:04random_translation_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype027
5random_translation_1/stateful_uniform/StatefulUniform�
)random_translation_1/stateful_uniform/subSub2random_translation_1/stateful_uniform/max:output:02random_translation_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2+
)random_translation_1/stateful_uniform/sub�
)random_translation_1/stateful_uniform/mulMul>random_translation_1/stateful_uniform/StatefulUniform:output:0-random_translation_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������2+
)random_translation_1/stateful_uniform/mul�
%random_translation_1/stateful_uniformAdd-random_translation_1/stateful_uniform/mul:z:02random_translation_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2'
%random_translation_1/stateful_uniform�
random_translation_1/mulMul)random_translation_1/stateful_uniform:z:0random_translation_1/Cast:y:0*
T0*'
_output_shapes
:���������2
random_translation_1/mul�
/random_translation_1/stateful_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :21
/random_translation_1/stateful_uniform_1/shape/1�
-random_translation_1/stateful_uniform_1/shapePack+random_translation_1/strided_slice:output:08random_translation_1/stateful_uniform_1/shape/1:output:0*
N*
T0*
_output_shapes
:2/
-random_translation_1/stateful_uniform_1/shape�
+random_translation_1/stateful_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2-
+random_translation_1/stateful_uniform_1/min�
+random_translation_1/stateful_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *    2-
+random_translation_1/stateful_uniform_1/max�
Arandom_translation_1/stateful_uniform_1/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2C
Arandom_translation_1/stateful_uniform_1/StatefulUniform/algorithm�
7random_translation_1/stateful_uniform_1/StatefulUniformStatefulUniform>random_translation_1_stateful_uniform_statefuluniform_resourceJrandom_translation_1/stateful_uniform_1/StatefulUniform/algorithm:output:06random_translation_1/stateful_uniform_1/shape:output:06^random_translation_1/stateful_uniform/StatefulUniform*'
_output_shapes
:���������*
shape_dtype029
7random_translation_1/stateful_uniform_1/StatefulUniform�
+random_translation_1/stateful_uniform_1/subSub4random_translation_1/stateful_uniform_1/max:output:04random_translation_1/stateful_uniform_1/min:output:0*
T0*
_output_shapes
: 2-
+random_translation_1/stateful_uniform_1/sub�
+random_translation_1/stateful_uniform_1/mulMul@random_translation_1/stateful_uniform_1/StatefulUniform:output:0/random_translation_1/stateful_uniform_1/sub:z:0*
T0*'
_output_shapes
:���������2-
+random_translation_1/stateful_uniform_1/mul�
'random_translation_1/stateful_uniform_1Add/random_translation_1/stateful_uniform_1/mul:z:04random_translation_1/stateful_uniform_1/min:output:0*
T0*'
_output_shapes
:���������2)
'random_translation_1/stateful_uniform_1�
random_translation_1/mul_1Mul+random_translation_1/stateful_uniform_1:z:0random_translation_1/Cast_1:y:0*
T0*'
_output_shapes
:���������2
random_translation_1/mul_1�
 random_translation_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2"
 random_translation_1/concat/axis�
random_translation_1/concatConcatV2random_translation_1/mul_1:z:0random_translation_1/mul:z:0)random_translation_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
random_translation_1/concat�
-random_translation_1/translation_matrix/ShapeShape$random_translation_1/concat:output:0*
T0*
_output_shapes
:2/
-random_translation_1/translation_matrix/Shape�
;random_translation_1/translation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2=
;random_translation_1/translation_matrix/strided_slice/stack�
=random_translation_1/translation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2?
=random_translation_1/translation_matrix/strided_slice/stack_1�
=random_translation_1/translation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2?
=random_translation_1/translation_matrix/strided_slice/stack_2�
5random_translation_1/translation_matrix/strided_sliceStridedSlice6random_translation_1/translation_matrix/Shape:output:0Drandom_translation_1/translation_matrix/strided_slice/stack:output:0Frandom_translation_1/translation_matrix/strided_slice/stack_1:output:0Frandom_translation_1/translation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask27
5random_translation_1/translation_matrix/strided_slice�
2random_translation_1/translation_matrix/ones/mul/yConst*
_output_shapes
: *
dtype0*
value	B :24
2random_translation_1/translation_matrix/ones/mul/y�
0random_translation_1/translation_matrix/ones/mulMul>random_translation_1/translation_matrix/strided_slice:output:0;random_translation_1/translation_matrix/ones/mul/y:output:0*
T0*
_output_shapes
: 22
0random_translation_1/translation_matrix/ones/mul�
3random_translation_1/translation_matrix/ones/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�25
3random_translation_1/translation_matrix/ones/Less/y�
1random_translation_1/translation_matrix/ones/LessLess4random_translation_1/translation_matrix/ones/mul:z:0<random_translation_1/translation_matrix/ones/Less/y:output:0*
T0*
_output_shapes
: 23
1random_translation_1/translation_matrix/ones/Less�
5random_translation_1/translation_matrix/ones/packed/1Const*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/ones/packed/1�
3random_translation_1/translation_matrix/ones/packedPack>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/ones/packed/1:output:0*
N*
T0*
_output_shapes
:25
3random_translation_1/translation_matrix/ones/packed�
2random_translation_1/translation_matrix/ones/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?24
2random_translation_1/translation_matrix/ones/Const�
,random_translation_1/translation_matrix/onesFill<random_translation_1/translation_matrix/ones/packed:output:0;random_translation_1/translation_matrix/ones/Const:output:0*
T0*'
_output_shapes
:���������2.
,random_translation_1/translation_matrix/ones�
3random_translation_1/translation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :25
3random_translation_1/translation_matrix/zeros/mul/y�
1random_translation_1/translation_matrix/zeros/mulMul>random_translation_1/translation_matrix/strided_slice:output:0<random_translation_1/translation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 23
1random_translation_1/translation_matrix/zeros/mul�
4random_translation_1/translation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�26
4random_translation_1/translation_matrix/zeros/Less/y�
2random_translation_1/translation_matrix/zeros/LessLess5random_translation_1/translation_matrix/zeros/mul:z:0=random_translation_1/translation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 24
2random_translation_1/translation_matrix/zeros/Less�
6random_translation_1/translation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :28
6random_translation_1/translation_matrix/zeros/packed/1�
4random_translation_1/translation_matrix/zeros/packedPack>random_translation_1/translation_matrix/strided_slice:output:0?random_translation_1/translation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:26
4random_translation_1/translation_matrix/zeros/packed�
3random_translation_1/translation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    25
3random_translation_1/translation_matrix/zeros/Const�
-random_translation_1/translation_matrix/zerosFill=random_translation_1/translation_matrix/zeros/packed:output:0<random_translation_1/translation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2/
-random_translation_1/translation_matrix/zeros�
=random_translation_1/translation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2?
=random_translation_1/translation_matrix/strided_slice_1/stack�
?random_translation_1/translation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2A
?random_translation_1/translation_matrix/strided_slice_1/stack_1�
?random_translation_1/translation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2A
?random_translation_1/translation_matrix/strided_slice_1/stack_2�
7random_translation_1/translation_matrix/strided_slice_1StridedSlice$random_translation_1/concat:output:0Frandom_translation_1/translation_matrix/strided_slice_1/stack:output:0Hrandom_translation_1/translation_matrix/strided_slice_1/stack_1:output:0Hrandom_translation_1/translation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask29
7random_translation_1/translation_matrix/strided_slice_1�
+random_translation_1/translation_matrix/NegNeg@random_translation_1/translation_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2-
+random_translation_1/translation_matrix/Neg�
5random_translation_1/translation_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/zeros_1/mul/y�
3random_translation_1/translation_matrix/zeros_1/mulMul>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/zeros_1/mul�
6random_translation_1/translation_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�28
6random_translation_1/translation_matrix/zeros_1/Less/y�
4random_translation_1/translation_matrix/zeros_1/LessLess7random_translation_1/translation_matrix/zeros_1/mul:z:0?random_translation_1/translation_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 26
4random_translation_1/translation_matrix/zeros_1/Less�
8random_translation_1/translation_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2:
8random_translation_1/translation_matrix/zeros_1/packed/1�
6random_translation_1/translation_matrix/zeros_1/packedPack>random_translation_1/translation_matrix/strided_slice:output:0Arandom_translation_1/translation_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:28
6random_translation_1/translation_matrix/zeros_1/packed�
5random_translation_1/translation_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    27
5random_translation_1/translation_matrix/zeros_1/Const�
/random_translation_1/translation_matrix/zeros_1Fill?random_translation_1/translation_matrix/zeros_1/packed:output:0>random_translation_1/translation_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������21
/random_translation_1/translation_matrix/zeros_1�
4random_translation_1/translation_matrix/ones_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :26
4random_translation_1/translation_matrix/ones_1/mul/y�
2random_translation_1/translation_matrix/ones_1/mulMul>random_translation_1/translation_matrix/strided_slice:output:0=random_translation_1/translation_matrix/ones_1/mul/y:output:0*
T0*
_output_shapes
: 24
2random_translation_1/translation_matrix/ones_1/mul�
5random_translation_1/translation_matrix/ones_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�27
5random_translation_1/translation_matrix/ones_1/Less/y�
3random_translation_1/translation_matrix/ones_1/LessLess6random_translation_1/translation_matrix/ones_1/mul:z:0>random_translation_1/translation_matrix/ones_1/Less/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/ones_1/Less�
7random_translation_1/translation_matrix/ones_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :29
7random_translation_1/translation_matrix/ones_1/packed/1�
5random_translation_1/translation_matrix/ones_1/packedPack>random_translation_1/translation_matrix/strided_slice:output:0@random_translation_1/translation_matrix/ones_1/packed/1:output:0*
N*
T0*
_output_shapes
:27
5random_translation_1/translation_matrix/ones_1/packed�
4random_translation_1/translation_matrix/ones_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  �?26
4random_translation_1/translation_matrix/ones_1/Const�
.random_translation_1/translation_matrix/ones_1Fill>random_translation_1/translation_matrix/ones_1/packed:output:0=random_translation_1/translation_matrix/ones_1/Const:output:0*
T0*'
_output_shapes
:���������20
.random_translation_1/translation_matrix/ones_1�
=random_translation_1/translation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2?
=random_translation_1/translation_matrix/strided_slice_2/stack�
?random_translation_1/translation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2A
?random_translation_1/translation_matrix/strided_slice_2/stack_1�
?random_translation_1/translation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2A
?random_translation_1/translation_matrix/strided_slice_2/stack_2�
7random_translation_1/translation_matrix/strided_slice_2StridedSlice$random_translation_1/concat:output:0Frandom_translation_1/translation_matrix/strided_slice_2/stack:output:0Hrandom_translation_1/translation_matrix/strided_slice_2/stack_1:output:0Hrandom_translation_1/translation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask29
7random_translation_1/translation_matrix/strided_slice_2�
-random_translation_1/translation_matrix/Neg_1Neg@random_translation_1/translation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2/
-random_translation_1/translation_matrix/Neg_1�
5random_translation_1/translation_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :27
5random_translation_1/translation_matrix/zeros_2/mul/y�
3random_translation_1/translation_matrix/zeros_2/mulMul>random_translation_1/translation_matrix/strided_slice:output:0>random_translation_1/translation_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 25
3random_translation_1/translation_matrix/zeros_2/mul�
6random_translation_1/translation_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�28
6random_translation_1/translation_matrix/zeros_2/Less/y�
4random_translation_1/translation_matrix/zeros_2/LessLess7random_translation_1/translation_matrix/zeros_2/mul:z:0?random_translation_1/translation_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 26
4random_translation_1/translation_matrix/zeros_2/Less�
8random_translation_1/translation_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2:
8random_translation_1/translation_matrix/zeros_2/packed/1�
6random_translation_1/translation_matrix/zeros_2/packedPack>random_translation_1/translation_matrix/strided_slice:output:0Arandom_translation_1/translation_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:28
6random_translation_1/translation_matrix/zeros_2/packed�
5random_translation_1/translation_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    27
5random_translation_1/translation_matrix/zeros_2/Const�
/random_translation_1/translation_matrix/zeros_2Fill?random_translation_1/translation_matrix/zeros_2/packed:output:0>random_translation_1/translation_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������21
/random_translation_1/translation_matrix/zeros_2�
3random_translation_1/translation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :25
3random_translation_1/translation_matrix/concat/axis�
.random_translation_1/translation_matrix/concatConcatV25random_translation_1/translation_matrix/ones:output:06random_translation_1/translation_matrix/zeros:output:0/random_translation_1/translation_matrix/Neg:y:08random_translation_1/translation_matrix/zeros_1:output:07random_translation_1/translation_matrix/ones_1:output:01random_translation_1/translation_matrix/Neg_1:y:08random_translation_1/translation_matrix/zeros_2:output:0<random_translation_1/translation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������20
.random_translation_1/translation_matrix/concat�
$random_translation_1/transform/ShapeShapeKrandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2&
$random_translation_1/transform/Shape�
2random_translation_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:24
2random_translation_1/transform/strided_slice/stack�
4random_translation_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:26
4random_translation_1/transform/strided_slice/stack_1�
4random_translation_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:26
4random_translation_1/transform/strided_slice/stack_2�
,random_translation_1/transform/strided_sliceStridedSlice-random_translation_1/transform/Shape:output:0;random_translation_1/transform/strided_slice/stack:output:0=random_translation_1/transform/strided_slice/stack_1:output:0=random_translation_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2.
,random_translation_1/transform/strided_slice�
)random_translation_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2+
)random_translation_1/transform/fill_value�
9random_translation_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Krandom_rotation_1/transform/ImageProjectiveTransformV3:transformed_images:07random_translation_1/translation_matrix/concat:output:05random_translation_1/transform/strided_slice:output:02random_translation_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	NEAREST*
interpolation
BILINEAR2;
9random_translation_1/transform/ImageProjectiveTransformV3�
random_zoom_1/ShapeShapeNrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_zoom_1/Shape�
!random_zoom_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2#
!random_zoom_1/strided_slice/stack�
#random_zoom_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice/stack_1�
#random_zoom_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice/stack_2�
random_zoom_1/strided_sliceStridedSlicerandom_zoom_1/Shape:output:0*random_zoom_1/strided_slice/stack:output:0,random_zoom_1/strided_slice/stack_1:output:0,random_zoom_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice�
#random_zoom_1/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice_1/stack�
%random_zoom_1/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_1/stack_1�
%random_zoom_1/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_1/stack_2�
random_zoom_1/strided_slice_1StridedSlicerandom_zoom_1/Shape:output:0,random_zoom_1/strided_slice_1/stack:output:0.random_zoom_1/strided_slice_1/stack_1:output:0.random_zoom_1/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice_1�
random_zoom_1/CastCast&random_zoom_1/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_zoom_1/Cast�
#random_zoom_1/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom_1/strided_slice_2/stack�
%random_zoom_1/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_2/stack_1�
%random_zoom_1/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_zoom_1/strided_slice_2/stack_2�
random_zoom_1/strided_slice_2StridedSlicerandom_zoom_1/Shape:output:0,random_zoom_1/strided_slice_2/stack:output:0.random_zoom_1/strided_slice_2/stack_1:output:0.random_zoom_1/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom_1/strided_slice_2�
random_zoom_1/Cast_1Cast&random_zoom_1/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_zoom_1/Cast_1�
&random_zoom_1/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2(
&random_zoom_1/stateful_uniform/shape/1�
$random_zoom_1/stateful_uniform/shapePack$random_zoom_1/strided_slice:output:0/random_zoom_1/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2&
$random_zoom_1/stateful_uniform/shape�
"random_zoom_1/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *��L?2$
"random_zoom_1/stateful_uniform/min�
"random_zoom_1/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���?2$
"random_zoom_1/stateful_uniform/max�
8random_zoom_1/stateful_uniform/StatefulUniform/algorithmConst*
_output_shapes
: *
dtype0	*
value	B	 R2:
8random_zoom_1/stateful_uniform/StatefulUniform/algorithm�
.random_zoom_1/stateful_uniform/StatefulUniformStatefulUniform7random_zoom_1_stateful_uniform_statefuluniform_resourceArandom_zoom_1/stateful_uniform/StatefulUniform/algorithm:output:0-random_zoom_1/stateful_uniform/shape:output:0*'
_output_shapes
:���������*
shape_dtype020
.random_zoom_1/stateful_uniform/StatefulUniform�
"random_zoom_1/stateful_uniform/subSub+random_zoom_1/stateful_uniform/max:output:0+random_zoom_1/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2$
"random_zoom_1/stateful_uniform/sub�
"random_zoom_1/stateful_uniform/mulMul7random_zoom_1/stateful_uniform/StatefulUniform:output:0&random_zoom_1/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:���������2$
"random_zoom_1/stateful_uniform/mul�
random_zoom_1/stateful_uniformAdd&random_zoom_1/stateful_uniform/mul:z:0+random_zoom_1/stateful_uniform/min:output:0*
T0*'
_output_shapes
:���������2 
random_zoom_1/stateful_uniformx
random_zoom_1/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
random_zoom_1/concat/axis�
random_zoom_1/concatConcatV2"random_zoom_1/stateful_uniform:z:0"random_zoom_1/stateful_uniform:z:0"random_zoom_1/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2
random_zoom_1/concat�
random_zoom_1/zoom_matrix/ShapeShaperandom_zoom_1/concat:output:0*
T0*
_output_shapes
:2!
random_zoom_1/zoom_matrix/Shape�
-random_zoom_1/zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2/
-random_zoom_1/zoom_matrix/strided_slice/stack�
/random_zoom_1/zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:21
/random_zoom_1/zoom_matrix/strided_slice/stack_1�
/random_zoom_1/zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:21
/random_zoom_1/zoom_matrix/strided_slice/stack_2�
'random_zoom_1/zoom_matrix/strided_sliceStridedSlice(random_zoom_1/zoom_matrix/Shape:output:06random_zoom_1/zoom_matrix/strided_slice/stack:output:08random_zoom_1/zoom_matrix/strided_slice/stack_1:output:08random_zoom_1/zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2)
'random_zoom_1/zoom_matrix/strided_slice�
random_zoom_1/zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2!
random_zoom_1/zoom_matrix/sub/y�
random_zoom_1/zoom_matrix/subSubrandom_zoom_1/Cast_1:y:0(random_zoom_1/zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 2
random_zoom_1/zoom_matrix/sub�
#random_zoom_1/zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2%
#random_zoom_1/zoom_matrix/truediv/y�
!random_zoom_1/zoom_matrix/truedivRealDiv!random_zoom_1/zoom_matrix/sub:z:0,random_zoom_1/zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 2#
!random_zoom_1/zoom_matrix/truediv�
/random_zoom_1/zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            21
/random_zoom_1/zoom_matrix/strided_slice_1/stack�
1random_zoom_1/zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_1/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_1/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_1StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_1/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_1/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_1�
!random_zoom_1/zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_1/x�
random_zoom_1/zoom_matrix/sub_1Sub*random_zoom_1/zoom_matrix/sub_1/x:output:02random_zoom_1/zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/sub_1�
random_zoom_1/zoom_matrix/mulMul%random_zoom_1/zoom_matrix/truediv:z:0#random_zoom_1/zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:���������2
random_zoom_1/zoom_matrix/mul�
!random_zoom_1/zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_2/y�
random_zoom_1/zoom_matrix/sub_2Subrandom_zoom_1/Cast:y:0*random_zoom_1/zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2!
random_zoom_1/zoom_matrix/sub_2�
%random_zoom_1/zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2'
%random_zoom_1/zoom_matrix/truediv_1/y�
#random_zoom_1/zoom_matrix/truediv_1RealDiv#random_zoom_1/zoom_matrix/sub_2:z:0.random_zoom_1/zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 2%
#random_zoom_1/zoom_matrix/truediv_1�
/random_zoom_1/zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom_1/zoom_matrix/strided_slice_2/stack�
1random_zoom_1/zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_2/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_2/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_2StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_2/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_2/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_2�
!random_zoom_1/zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?2#
!random_zoom_1/zoom_matrix/sub_3/x�
random_zoom_1/zoom_matrix/sub_3Sub*random_zoom_1/zoom_matrix/sub_3/x:output:02random_zoom_1/zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/sub_3�
random_zoom_1/zoom_matrix/mul_1Mul'random_zoom_1/zoom_matrix/truediv_1:z:0#random_zoom_1/zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/mul_1�
/random_zoom_1/zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            21
/random_zoom_1/zoom_matrix/strided_slice_3/stack�
1random_zoom_1/zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_3/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_3/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_3StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_3/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_3/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_3�
%random_zoom_1/zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2'
%random_zoom_1/zoom_matrix/zeros/mul/y�
#random_zoom_1/zoom_matrix/zeros/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:0.random_zoom_1/zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2%
#random_zoom_1/zoom_matrix/zeros/mul�
&random_zoom_1/zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2(
&random_zoom_1/zoom_matrix/zeros/Less/y�
$random_zoom_1/zoom_matrix/zeros/LessLess'random_zoom_1/zoom_matrix/zeros/mul:z:0/random_zoom_1/zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2&
$random_zoom_1/zoom_matrix/zeros/Less�
(random_zoom_1/zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2*
(random_zoom_1/zoom_matrix/zeros/packed/1�
&random_zoom_1/zoom_matrix/zeros/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:01random_zoom_1/zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2(
&random_zoom_1/zoom_matrix/zeros/packed�
%random_zoom_1/zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%random_zoom_1/zoom_matrix/zeros/Const�
random_zoom_1/zoom_matrix/zerosFill/random_zoom_1/zoom_matrix/zeros/packed:output:0.random_zoom_1/zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:���������2!
random_zoom_1/zoom_matrix/zeros�
'random_zoom_1/zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2)
'random_zoom_1/zoom_matrix/zeros_1/mul/y�
%random_zoom_1/zoom_matrix/zeros_1/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:00random_zoom_1/zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2'
%random_zoom_1/zoom_matrix/zeros_1/mul�
(random_zoom_1/zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2*
(random_zoom_1/zoom_matrix/zeros_1/Less/y�
&random_zoom_1/zoom_matrix/zeros_1/LessLess)random_zoom_1/zoom_matrix/zeros_1/mul:z:01random_zoom_1/zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2(
&random_zoom_1/zoom_matrix/zeros_1/Less�
*random_zoom_1/zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2,
*random_zoom_1/zoom_matrix/zeros_1/packed/1�
(random_zoom_1/zoom_matrix/zeros_1/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:03random_zoom_1/zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2*
(random_zoom_1/zoom_matrix/zeros_1/packed�
'random_zoom_1/zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2)
'random_zoom_1/zoom_matrix/zeros_1/Const�
!random_zoom_1/zoom_matrix/zeros_1Fill1random_zoom_1/zoom_matrix/zeros_1/packed:output:00random_zoom_1/zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:���������2#
!random_zoom_1/zoom_matrix/zeros_1�
/random_zoom_1/zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom_1/zoom_matrix/strided_slice_4/stack�
1random_zoom_1/zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           23
1random_zoom_1/zoom_matrix/strided_slice_4/stack_1�
1random_zoom_1/zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         23
1random_zoom_1/zoom_matrix/strided_slice_4/stack_2�
)random_zoom_1/zoom_matrix/strided_slice_4StridedSlicerandom_zoom_1/concat:output:08random_zoom_1/zoom_matrix/strided_slice_4/stack:output:0:random_zoom_1/zoom_matrix/strided_slice_4/stack_1:output:0:random_zoom_1/zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:���������*

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2+
)random_zoom_1/zoom_matrix/strided_slice_4�
'random_zoom_1/zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2)
'random_zoom_1/zoom_matrix/zeros_2/mul/y�
%random_zoom_1/zoom_matrix/zeros_2/mulMul0random_zoom_1/zoom_matrix/strided_slice:output:00random_zoom_1/zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2'
%random_zoom_1/zoom_matrix/zeros_2/mul�
(random_zoom_1/zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :�2*
(random_zoom_1/zoom_matrix/zeros_2/Less/y�
&random_zoom_1/zoom_matrix/zeros_2/LessLess)random_zoom_1/zoom_matrix/zeros_2/mul:z:01random_zoom_1/zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2(
&random_zoom_1/zoom_matrix/zeros_2/Less�
*random_zoom_1/zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2,
*random_zoom_1/zoom_matrix/zeros_2/packed/1�
(random_zoom_1/zoom_matrix/zeros_2/packedPack0random_zoom_1/zoom_matrix/strided_slice:output:03random_zoom_1/zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2*
(random_zoom_1/zoom_matrix/zeros_2/packed�
'random_zoom_1/zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2)
'random_zoom_1/zoom_matrix/zeros_2/Const�
!random_zoom_1/zoom_matrix/zeros_2Fill1random_zoom_1/zoom_matrix/zeros_2/packed:output:00random_zoom_1/zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:���������2#
!random_zoom_1/zoom_matrix/zeros_2�
%random_zoom_1/zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2'
%random_zoom_1/zoom_matrix/concat/axis�
 random_zoom_1/zoom_matrix/concatConcatV22random_zoom_1/zoom_matrix/strided_slice_3:output:0(random_zoom_1/zoom_matrix/zeros:output:0!random_zoom_1/zoom_matrix/mul:z:0*random_zoom_1/zoom_matrix/zeros_1:output:02random_zoom_1/zoom_matrix/strided_slice_4:output:0#random_zoom_1/zoom_matrix/mul_1:z:0*random_zoom_1/zoom_matrix/zeros_2:output:0.random_zoom_1/zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:���������2"
 random_zoom_1/zoom_matrix/concat�
random_zoom_1/transform/ShapeShapeNrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_zoom_1/transform/Shape�
+random_zoom_1/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2-
+random_zoom_1/transform/strided_slice/stack�
-random_zoom_1/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2/
-random_zoom_1/transform/strided_slice/stack_1�
-random_zoom_1/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2/
-random_zoom_1/transform/strided_slice/stack_2�
%random_zoom_1/transform/strided_sliceStridedSlice&random_zoom_1/transform/Shape:output:04random_zoom_1/transform/strided_slice/stack:output:06random_zoom_1/transform/strided_slice/stack_1:output:06random_zoom_1/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2'
%random_zoom_1/transform/strided_slice�
"random_zoom_1/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2$
"random_zoom_1/transform/fill_value�
2random_zoom_1/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Nrandom_translation_1/transform/ImageProjectiveTransformV3:transformed_images:0)random_zoom_1/zoom_matrix/concat:output:0.random_zoom_1/transform/strided_slice:output:0+random_zoom_1/transform/fill_value:output:0*/
_output_shapes
:���������&*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR24
2random_zoom_1/transform/ImageProjectiveTransformV3�
IdentityIdentityGrandom_zoom_1/transform/ImageProjectiveTransformV3:transformed_images:03^random_rotation_1/stateful_uniform/StatefulUniform6^random_translation_1/stateful_uniform/StatefulUniform8^random_translation_1/stateful_uniform_1/StatefulUniform/^random_zoom_1/stateful_uniform/StatefulUniform*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*:
_input_shapes)
':���������&:::2h
2random_rotation_1/stateful_uniform/StatefulUniform2random_rotation_1/stateful_uniform/StatefulUniform2n
5random_translation_1/stateful_uniform/StatefulUniform5random_translation_1/stateful_uniform/StatefulUniform2r
7random_translation_1/stateful_uniform_1/StatefulUniform7random_translation_1/stateful_uniform_1/StatefulUniform2`
.random_zoom_1/stateful_uniform/StatefulUniform.random_zoom_1/stateful_uniform/StatefulUniform:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�s
�
 __inference__wrapped_model_30629
sequential_4_inputE
Asequential_5_sequential_4_conv2d_3_conv2d_readvariableop_resourceF
Bsequential_5_sequential_4_conv2d_3_biasadd_readvariableop_resourceE
Asequential_5_sequential_4_conv2d_4_conv2d_readvariableop_resourceF
Bsequential_5_sequential_4_conv2d_4_biasadd_readvariableop_resourceE
Asequential_5_sequential_4_conv2d_5_conv2d_readvariableop_resourceF
Bsequential_5_sequential_4_conv2d_5_biasadd_readvariableop_resourceD
@sequential_5_sequential_4_dense_3_matmul_readvariableop_resourceE
Asequential_5_sequential_4_dense_3_biasadd_readvariableop_resourceD
@sequential_5_sequential_4_dense_4_matmul_readvariableop_resourceE
Asequential_5_sequential_4_dense_4_biasadd_readvariableop_resourceD
@sequential_5_sequential_4_dense_5_matmul_readvariableop_resourceE
Asequential_5_sequential_4_dense_5_biasadd_readvariableop_resource
identity��9sequential_5/sequential_4/conv2d_3/BiasAdd/ReadVariableOp�8sequential_5/sequential_4/conv2d_3/Conv2D/ReadVariableOp�9sequential_5/sequential_4/conv2d_4/BiasAdd/ReadVariableOp�8sequential_5/sequential_4/conv2d_4/Conv2D/ReadVariableOp�9sequential_5/sequential_4/conv2d_5/BiasAdd/ReadVariableOp�8sequential_5/sequential_4/conv2d_5/Conv2D/ReadVariableOp�8sequential_5/sequential_4/dense_3/BiasAdd/ReadVariableOp�7sequential_5/sequential_4/dense_3/MatMul/ReadVariableOp�8sequential_5/sequential_4/dense_4/BiasAdd/ReadVariableOp�7sequential_5/sequential_4/dense_4/MatMul/ReadVariableOp�8sequential_5/sequential_4/dense_5/BiasAdd/ReadVariableOp�7sequential_5/sequential_4/dense_5/MatMul/ReadVariableOp�
,sequential_5/sequential_4/rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2.
,sequential_5/sequential_4/rescaling_3/Cast/x�
.sequential_5/sequential_4/rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    20
.sequential_5/sequential_4/rescaling_3/Cast_1/x�
)sequential_5/sequential_4/rescaling_3/mulMulsequential_4_input5sequential_5/sequential_4/rescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2+
)sequential_5/sequential_4/rescaling_3/mul�
)sequential_5/sequential_4/rescaling_3/addAddV2-sequential_5/sequential_4/rescaling_3/mul:z:07sequential_5/sequential_4/rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2+
)sequential_5/sequential_4/rescaling_3/add�
8sequential_5/sequential_4/conv2d_3/Conv2D/ReadVariableOpReadVariableOpAsequential_5_sequential_4_conv2d_3_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02:
8sequential_5/sequential_4/conv2d_3/Conv2D/ReadVariableOp�
)sequential_5/sequential_4/conv2d_3/Conv2DConv2D-sequential_5/sequential_4/rescaling_3/add:z:0@sequential_5/sequential_4/conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&*
paddingSAME*
strides
2+
)sequential_5/sequential_4/conv2d_3/Conv2D�
9sequential_5/sequential_4/conv2d_3/BiasAdd/ReadVariableOpReadVariableOpBsequential_5_sequential_4_conv2d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9sequential_5/sequential_4/conv2d_3/BiasAdd/ReadVariableOp�
*sequential_5/sequential_4/conv2d_3/BiasAddBiasAdd2sequential_5/sequential_4/conv2d_3/Conv2D:output:0Asequential_5/sequential_4/conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&2,
*sequential_5/sequential_4/conv2d_3/BiasAdd�
'sequential_5/sequential_4/conv2d_3/ReluRelu3sequential_5/sequential_4/conv2d_3/BiasAdd:output:0*
T0*/
_output_shapes
:���������&2)
'sequential_5/sequential_4/conv2d_3/Relu�
,sequential_5/sequential_4/dropout_4/IdentityIdentity5sequential_5/sequential_4/conv2d_3/Relu:activations:0*
T0*/
_output_shapes
:���������&2.
,sequential_5/sequential_4/dropout_4/Identity�
1sequential_5/sequential_4/max_pooling2d_3/MaxPoolMaxPool5sequential_5/sequential_4/dropout_4/Identity:output:0*/
_output_shapes
:���������*
ksize
*
paddingVALID*
strides
23
1sequential_5/sequential_4/max_pooling2d_3/MaxPool�
8sequential_5/sequential_4/conv2d_4/Conv2D/ReadVariableOpReadVariableOpAsequential_5_sequential_4_conv2d_4_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02:
8sequential_5/sequential_4/conv2d_4/Conv2D/ReadVariableOp�
)sequential_5/sequential_4/conv2d_4/Conv2DConv2D:sequential_5/sequential_4/max_pooling2d_3/MaxPool:output:0@sequential_5/sequential_4/conv2d_4/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������*
paddingSAME*
strides
2+
)sequential_5/sequential_4/conv2d_4/Conv2D�
9sequential_5/sequential_4/conv2d_4/BiasAdd/ReadVariableOpReadVariableOpBsequential_5_sequential_4_conv2d_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9sequential_5/sequential_4/conv2d_4/BiasAdd/ReadVariableOp�
*sequential_5/sequential_4/conv2d_4/BiasAddBiasAdd2sequential_5/sequential_4/conv2d_4/Conv2D:output:0Asequential_5/sequential_4/conv2d_4/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������2,
*sequential_5/sequential_4/conv2d_4/BiasAdd�
'sequential_5/sequential_4/conv2d_4/ReluRelu3sequential_5/sequential_4/conv2d_4/BiasAdd:output:0*
T0*/
_output_shapes
:���������2)
'sequential_5/sequential_4/conv2d_4/Relu�
,sequential_5/sequential_4/dropout_5/IdentityIdentity5sequential_5/sequential_4/conv2d_4/Relu:activations:0*
T0*/
_output_shapes
:���������2.
,sequential_5/sequential_4/dropout_5/Identity�
1sequential_5/sequential_4/max_pooling2d_4/MaxPoolMaxPool5sequential_5/sequential_4/dropout_5/Identity:output:0*/
_output_shapes
:���������	*
ksize
*
paddingVALID*
strides
23
1sequential_5/sequential_4/max_pooling2d_4/MaxPool�
8sequential_5/sequential_4/conv2d_5/Conv2D/ReadVariableOpReadVariableOpAsequential_5_sequential_4_conv2d_5_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype02:
8sequential_5/sequential_4/conv2d_5/Conv2D/ReadVariableOp�
)sequential_5/sequential_4/conv2d_5/Conv2DConv2D:sequential_5/sequential_4/max_pooling2d_4/MaxPool:output:0@sequential_5/sequential_4/conv2d_5/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 *
paddingSAME*
strides
2+
)sequential_5/sequential_4/conv2d_5/Conv2D�
9sequential_5/sequential_4/conv2d_5/BiasAdd/ReadVariableOpReadVariableOpBsequential_5_sequential_4_conv2d_5_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02;
9sequential_5/sequential_4/conv2d_5/BiasAdd/ReadVariableOp�
*sequential_5/sequential_4/conv2d_5/BiasAddBiasAdd2sequential_5/sequential_4/conv2d_5/Conv2D:output:0Asequential_5/sequential_4/conv2d_5/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 2,
*sequential_5/sequential_4/conv2d_5/BiasAdd�
'sequential_5/sequential_4/conv2d_5/ReluRelu3sequential_5/sequential_4/conv2d_5/BiasAdd:output:0*
T0*/
_output_shapes
:���������	 2)
'sequential_5/sequential_4/conv2d_5/Relu�
1sequential_5/sequential_4/max_pooling2d_5/MaxPoolMaxPool5sequential_5/sequential_4/conv2d_5/Relu:activations:0*/
_output_shapes
:��������� *
ksize
*
paddingVALID*
strides
23
1sequential_5/sequential_4/max_pooling2d_5/MaxPool�
)sequential_5/sequential_4/flatten_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"�����  2+
)sequential_5/sequential_4/flatten_1/Const�
+sequential_5/sequential_4/flatten_1/ReshapeReshape:sequential_5/sequential_4/max_pooling2d_5/MaxPool:output:02sequential_5/sequential_4/flatten_1/Const:output:0*
T0*(
_output_shapes
:����������2-
+sequential_5/sequential_4/flatten_1/Reshape�
7sequential_5/sequential_4/dense_3/MatMul/ReadVariableOpReadVariableOp@sequential_5_sequential_4_dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype029
7sequential_5/sequential_4/dense_3/MatMul/ReadVariableOp�
(sequential_5/sequential_4/dense_3/MatMulMatMul4sequential_5/sequential_4/flatten_1/Reshape:output:0?sequential_5/sequential_4/dense_3/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2*
(sequential_5/sequential_4/dense_3/MatMul�
8sequential_5/sequential_4/dense_3/BiasAdd/ReadVariableOpReadVariableOpAsequential_5_sequential_4_dense_3_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype02:
8sequential_5/sequential_4/dense_3/BiasAdd/ReadVariableOp�
)sequential_5/sequential_4/dense_3/BiasAddBiasAdd2sequential_5/sequential_4/dense_3/MatMul:product:0@sequential_5/sequential_4/dense_3/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2+
)sequential_5/sequential_4/dense_3/BiasAdd�
&sequential_5/sequential_4/dense_3/ReluRelu2sequential_5/sequential_4/dense_3/BiasAdd:output:0*
T0*(
_output_shapes
:����������2(
&sequential_5/sequential_4/dense_3/Relu�
,sequential_5/sequential_4/dropout_6/IdentityIdentity4sequential_5/sequential_4/dense_3/Relu:activations:0*
T0*(
_output_shapes
:����������2.
,sequential_5/sequential_4/dropout_6/Identity�
7sequential_5/sequential_4/dense_4/MatMul/ReadVariableOpReadVariableOp@sequential_5_sequential_4_dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype029
7sequential_5/sequential_4/dense_4/MatMul/ReadVariableOp�
(sequential_5/sequential_4/dense_4/MatMulMatMul5sequential_5/sequential_4/dropout_6/Identity:output:0?sequential_5/sequential_4/dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2*
(sequential_5/sequential_4/dense_4/MatMul�
8sequential_5/sequential_4/dense_4/BiasAdd/ReadVariableOpReadVariableOpAsequential_5_sequential_4_dense_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02:
8sequential_5/sequential_4/dense_4/BiasAdd/ReadVariableOp�
)sequential_5/sequential_4/dense_4/BiasAddBiasAdd2sequential_5/sequential_4/dense_4/MatMul:product:0@sequential_5/sequential_4/dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2+
)sequential_5/sequential_4/dense_4/BiasAdd�
&sequential_5/sequential_4/dense_4/ReluRelu2sequential_5/sequential_4/dense_4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2(
&sequential_5/sequential_4/dense_4/Relu�
,sequential_5/sequential_4/dropout_7/IdentityIdentity4sequential_5/sequential_4/dense_4/Relu:activations:0*
T0*'
_output_shapes
:���������2.
,sequential_5/sequential_4/dropout_7/Identity�
7sequential_5/sequential_4/dense_5/MatMul/ReadVariableOpReadVariableOp@sequential_5_sequential_4_dense_5_matmul_readvariableop_resource*
_output_shapes

:*
dtype029
7sequential_5/sequential_4/dense_5/MatMul/ReadVariableOp�
(sequential_5/sequential_4/dense_5/MatMulMatMul5sequential_5/sequential_4/dropout_7/Identity:output:0?sequential_5/sequential_4/dense_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2*
(sequential_5/sequential_4/dense_5/MatMul�
8sequential_5/sequential_4/dense_5/BiasAdd/ReadVariableOpReadVariableOpAsequential_5_sequential_4_dense_5_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02:
8sequential_5/sequential_4/dense_5/BiasAdd/ReadVariableOp�
)sequential_5/sequential_4/dense_5/BiasAddBiasAdd2sequential_5/sequential_4/dense_5/MatMul:product:0@sequential_5/sequential_4/dense_5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2+
)sequential_5/sequential_4/dense_5/BiasAdd�
)sequential_5/sequential_4/dense_5/SigmoidSigmoid2sequential_5/sequential_4/dense_5/BiasAdd:output:0*
T0*'
_output_shapes
:���������2+
)sequential_5/sequential_4/dense_5/Sigmoid�
sequential_5/softmax_1/SoftmaxSoftmax-sequential_5/sequential_4/dense_5/Sigmoid:y:0*
T0*'
_output_shapes
:���������2 
sequential_5/softmax_1/Softmax�
IdentityIdentity(sequential_5/softmax_1/Softmax:softmax:0:^sequential_5/sequential_4/conv2d_3/BiasAdd/ReadVariableOp9^sequential_5/sequential_4/conv2d_3/Conv2D/ReadVariableOp:^sequential_5/sequential_4/conv2d_4/BiasAdd/ReadVariableOp9^sequential_5/sequential_4/conv2d_4/Conv2D/ReadVariableOp:^sequential_5/sequential_4/conv2d_5/BiasAdd/ReadVariableOp9^sequential_5/sequential_4/conv2d_5/Conv2D/ReadVariableOp9^sequential_5/sequential_4/dense_3/BiasAdd/ReadVariableOp8^sequential_5/sequential_4/dense_3/MatMul/ReadVariableOp9^sequential_5/sequential_4/dense_4/BiasAdd/ReadVariableOp8^sequential_5/sequential_4/dense_4/MatMul/ReadVariableOp9^sequential_5/sequential_4/dense_5/BiasAdd/ReadVariableOp8^sequential_5/sequential_4/dense_5/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::2v
9sequential_5/sequential_4/conv2d_3/BiasAdd/ReadVariableOp9sequential_5/sequential_4/conv2d_3/BiasAdd/ReadVariableOp2t
8sequential_5/sequential_4/conv2d_3/Conv2D/ReadVariableOp8sequential_5/sequential_4/conv2d_3/Conv2D/ReadVariableOp2v
9sequential_5/sequential_4/conv2d_4/BiasAdd/ReadVariableOp9sequential_5/sequential_4/conv2d_4/BiasAdd/ReadVariableOp2t
8sequential_5/sequential_4/conv2d_4/Conv2D/ReadVariableOp8sequential_5/sequential_4/conv2d_4/Conv2D/ReadVariableOp2v
9sequential_5/sequential_4/conv2d_5/BiasAdd/ReadVariableOp9sequential_5/sequential_4/conv2d_5/BiasAdd/ReadVariableOp2t
8sequential_5/sequential_4/conv2d_5/Conv2D/ReadVariableOp8sequential_5/sequential_4/conv2d_5/Conv2D/ReadVariableOp2t
8sequential_5/sequential_4/dense_3/BiasAdd/ReadVariableOp8sequential_5/sequential_4/dense_3/BiasAdd/ReadVariableOp2r
7sequential_5/sequential_4/dense_3/MatMul/ReadVariableOp7sequential_5/sequential_4/dense_3/MatMul/ReadVariableOp2t
8sequential_5/sequential_4/dense_4/BiasAdd/ReadVariableOp8sequential_5/sequential_4/dense_4/BiasAdd/ReadVariableOp2r
7sequential_5/sequential_4/dense_4/MatMul/ReadVariableOp7sequential_5/sequential_4/dense_4/MatMul/ReadVariableOp2t
8sequential_5/sequential_4/dense_5/BiasAdd/ReadVariableOp8sequential_5/sequential_4/dense_5/BiasAdd/ReadVariableOp2r
7sequential_5/sequential_4/dense_5/MatMul/ReadVariableOp7sequential_5/sequential_4/dense_5/MatMul/ReadVariableOp:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_4_input
�s
�
G__inference_sequential_5_layer_call_and_return_conditional_losses_32811

inputs8
4sequential_4_conv2d_3_conv2d_readvariableop_resource9
5sequential_4_conv2d_3_biasadd_readvariableop_resource8
4sequential_4_conv2d_4_conv2d_readvariableop_resource9
5sequential_4_conv2d_4_biasadd_readvariableop_resource8
4sequential_4_conv2d_5_conv2d_readvariableop_resource9
5sequential_4_conv2d_5_biasadd_readvariableop_resource7
3sequential_4_dense_3_matmul_readvariableop_resource8
4sequential_4_dense_3_biasadd_readvariableop_resource7
3sequential_4_dense_4_matmul_readvariableop_resource8
4sequential_4_dense_4_biasadd_readvariableop_resource7
3sequential_4_dense_5_matmul_readvariableop_resource8
4sequential_4_dense_5_biasadd_readvariableop_resource
identity��0dense_3/kernel/Regularizer/Square/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�,sequential_4/conv2d_3/BiasAdd/ReadVariableOp�+sequential_4/conv2d_3/Conv2D/ReadVariableOp�,sequential_4/conv2d_4/BiasAdd/ReadVariableOp�+sequential_4/conv2d_4/Conv2D/ReadVariableOp�,sequential_4/conv2d_5/BiasAdd/ReadVariableOp�+sequential_4/conv2d_5/Conv2D/ReadVariableOp�+sequential_4/dense_3/BiasAdd/ReadVariableOp�*sequential_4/dense_3/MatMul/ReadVariableOp�+sequential_4/dense_4/BiasAdd/ReadVariableOp�*sequential_4/dense_4/MatMul/ReadVariableOp�+sequential_4/dense_5/BiasAdd/ReadVariableOp�*sequential_4/dense_5/MatMul/ReadVariableOp�
sequential_4/rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2!
sequential_4/rescaling_3/Cast/x�
!sequential_4/rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!sequential_4/rescaling_3/Cast_1/x�
sequential_4/rescaling_3/mulMulinputs(sequential_4/rescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2
sequential_4/rescaling_3/mul�
sequential_4/rescaling_3/addAddV2 sequential_4/rescaling_3/mul:z:0*sequential_4/rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2
sequential_4/rescaling_3/add�
+sequential_4/conv2d_3/Conv2D/ReadVariableOpReadVariableOp4sequential_4_conv2d_3_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02-
+sequential_4/conv2d_3/Conv2D/ReadVariableOp�
sequential_4/conv2d_3/Conv2DConv2D sequential_4/rescaling_3/add:z:03sequential_4/conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&*
paddingSAME*
strides
2
sequential_4/conv2d_3/Conv2D�
,sequential_4/conv2d_3/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_conv2d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,sequential_4/conv2d_3/BiasAdd/ReadVariableOp�
sequential_4/conv2d_3/BiasAddBiasAdd%sequential_4/conv2d_3/Conv2D:output:04sequential_4/conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&2
sequential_4/conv2d_3/BiasAdd�
sequential_4/conv2d_3/ReluRelu&sequential_4/conv2d_3/BiasAdd:output:0*
T0*/
_output_shapes
:���������&2
sequential_4/conv2d_3/Relu�
sequential_4/dropout_4/IdentityIdentity(sequential_4/conv2d_3/Relu:activations:0*
T0*/
_output_shapes
:���������&2!
sequential_4/dropout_4/Identity�
$sequential_4/max_pooling2d_3/MaxPoolMaxPool(sequential_4/dropout_4/Identity:output:0*/
_output_shapes
:���������*
ksize
*
paddingVALID*
strides
2&
$sequential_4/max_pooling2d_3/MaxPool�
+sequential_4/conv2d_4/Conv2D/ReadVariableOpReadVariableOp4sequential_4_conv2d_4_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02-
+sequential_4/conv2d_4/Conv2D/ReadVariableOp�
sequential_4/conv2d_4/Conv2DConv2D-sequential_4/max_pooling2d_3/MaxPool:output:03sequential_4/conv2d_4/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������*
paddingSAME*
strides
2
sequential_4/conv2d_4/Conv2D�
,sequential_4/conv2d_4/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_conv2d_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,sequential_4/conv2d_4/BiasAdd/ReadVariableOp�
sequential_4/conv2d_4/BiasAddBiasAdd%sequential_4/conv2d_4/Conv2D:output:04sequential_4/conv2d_4/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������2
sequential_4/conv2d_4/BiasAdd�
sequential_4/conv2d_4/ReluRelu&sequential_4/conv2d_4/BiasAdd:output:0*
T0*/
_output_shapes
:���������2
sequential_4/conv2d_4/Relu�
sequential_4/dropout_5/IdentityIdentity(sequential_4/conv2d_4/Relu:activations:0*
T0*/
_output_shapes
:���������2!
sequential_4/dropout_5/Identity�
$sequential_4/max_pooling2d_4/MaxPoolMaxPool(sequential_4/dropout_5/Identity:output:0*/
_output_shapes
:���������	*
ksize
*
paddingVALID*
strides
2&
$sequential_4/max_pooling2d_4/MaxPool�
+sequential_4/conv2d_5/Conv2D/ReadVariableOpReadVariableOp4sequential_4_conv2d_5_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype02-
+sequential_4/conv2d_5/Conv2D/ReadVariableOp�
sequential_4/conv2d_5/Conv2DConv2D-sequential_4/max_pooling2d_4/MaxPool:output:03sequential_4/conv2d_5/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 *
paddingSAME*
strides
2
sequential_4/conv2d_5/Conv2D�
,sequential_4/conv2d_5/BiasAdd/ReadVariableOpReadVariableOp5sequential_4_conv2d_5_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02.
,sequential_4/conv2d_5/BiasAdd/ReadVariableOp�
sequential_4/conv2d_5/BiasAddBiasAdd%sequential_4/conv2d_5/Conv2D:output:04sequential_4/conv2d_5/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 2
sequential_4/conv2d_5/BiasAdd�
sequential_4/conv2d_5/ReluRelu&sequential_4/conv2d_5/BiasAdd:output:0*
T0*/
_output_shapes
:���������	 2
sequential_4/conv2d_5/Relu�
$sequential_4/max_pooling2d_5/MaxPoolMaxPool(sequential_4/conv2d_5/Relu:activations:0*/
_output_shapes
:��������� *
ksize
*
paddingVALID*
strides
2&
$sequential_4/max_pooling2d_5/MaxPool�
sequential_4/flatten_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"�����  2
sequential_4/flatten_1/Const�
sequential_4/flatten_1/ReshapeReshape-sequential_4/max_pooling2d_5/MaxPool:output:0%sequential_4/flatten_1/Const:output:0*
T0*(
_output_shapes
:����������2 
sequential_4/flatten_1/Reshape�
*sequential_4/dense_3/MatMul/ReadVariableOpReadVariableOp3sequential_4_dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype02,
*sequential_4/dense_3/MatMul/ReadVariableOp�
sequential_4/dense_3/MatMulMatMul'sequential_4/flatten_1/Reshape:output:02sequential_4/dense_3/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
sequential_4/dense_3/MatMul�
+sequential_4/dense_3/BiasAdd/ReadVariableOpReadVariableOp4sequential_4_dense_3_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype02-
+sequential_4/dense_3/BiasAdd/ReadVariableOp�
sequential_4/dense_3/BiasAddBiasAdd%sequential_4/dense_3/MatMul:product:03sequential_4/dense_3/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
sequential_4/dense_3/BiasAdd�
sequential_4/dense_3/ReluRelu%sequential_4/dense_3/BiasAdd:output:0*
T0*(
_output_shapes
:����������2
sequential_4/dense_3/Relu�
sequential_4/dropout_6/IdentityIdentity'sequential_4/dense_3/Relu:activations:0*
T0*(
_output_shapes
:����������2!
sequential_4/dropout_6/Identity�
*sequential_4/dense_4/MatMul/ReadVariableOpReadVariableOp3sequential_4_dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype02,
*sequential_4/dense_4/MatMul/ReadVariableOp�
sequential_4/dense_4/MatMulMatMul(sequential_4/dropout_6/Identity:output:02sequential_4/dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_4/MatMul�
+sequential_4/dense_4/BiasAdd/ReadVariableOpReadVariableOp4sequential_4_dense_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+sequential_4/dense_4/BiasAdd/ReadVariableOp�
sequential_4/dense_4/BiasAddBiasAdd%sequential_4/dense_4/MatMul:product:03sequential_4/dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_4/BiasAdd�
sequential_4/dense_4/ReluRelu%sequential_4/dense_4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_4/Relu�
sequential_4/dropout_7/IdentityIdentity'sequential_4/dense_4/Relu:activations:0*
T0*'
_output_shapes
:���������2!
sequential_4/dropout_7/Identity�
*sequential_4/dense_5/MatMul/ReadVariableOpReadVariableOp3sequential_4_dense_5_matmul_readvariableop_resource*
_output_shapes

:*
dtype02,
*sequential_4/dense_5/MatMul/ReadVariableOp�
sequential_4/dense_5/MatMulMatMul(sequential_4/dropout_7/Identity:output:02sequential_4/dense_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_5/MatMul�
+sequential_4/dense_5/BiasAdd/ReadVariableOpReadVariableOp4sequential_4_dense_5_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+sequential_4/dense_5/BiasAdd/ReadVariableOp�
sequential_4/dense_5/BiasAddBiasAdd%sequential_4/dense_5/MatMul:product:03sequential_4/dense_5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_5/BiasAdd�
sequential_4/dense_5/SigmoidSigmoid%sequential_4/dense_5/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
sequential_4/dense_5/Sigmoid�
softmax_1/SoftmaxSoftmax sequential_4/dense_5/Sigmoid:y:0*
T0*'
_output_shapes
:���������2
softmax_1/Softmax�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp3sequential_4_dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOp3sequential_4_dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentitysoftmax_1/Softmax:softmax:01^dense_3/kernel/Regularizer/Square/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp-^sequential_4/conv2d_3/BiasAdd/ReadVariableOp,^sequential_4/conv2d_3/Conv2D/ReadVariableOp-^sequential_4/conv2d_4/BiasAdd/ReadVariableOp,^sequential_4/conv2d_4/Conv2D/ReadVariableOp-^sequential_4/conv2d_5/BiasAdd/ReadVariableOp,^sequential_4/conv2d_5/Conv2D/ReadVariableOp,^sequential_4/dense_3/BiasAdd/ReadVariableOp+^sequential_4/dense_3/MatMul/ReadVariableOp,^sequential_4/dense_4/BiasAdd/ReadVariableOp+^sequential_4/dense_4/MatMul/ReadVariableOp,^sequential_4/dense_5/BiasAdd/ReadVariableOp+^sequential_4/dense_5/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2\
,sequential_4/conv2d_3/BiasAdd/ReadVariableOp,sequential_4/conv2d_3/BiasAdd/ReadVariableOp2Z
+sequential_4/conv2d_3/Conv2D/ReadVariableOp+sequential_4/conv2d_3/Conv2D/ReadVariableOp2\
,sequential_4/conv2d_4/BiasAdd/ReadVariableOp,sequential_4/conv2d_4/BiasAdd/ReadVariableOp2Z
+sequential_4/conv2d_4/Conv2D/ReadVariableOp+sequential_4/conv2d_4/Conv2D/ReadVariableOp2\
,sequential_4/conv2d_5/BiasAdd/ReadVariableOp,sequential_4/conv2d_5/BiasAdd/ReadVariableOp2Z
+sequential_4/conv2d_5/Conv2D/ReadVariableOp+sequential_4/conv2d_5/Conv2D/ReadVariableOp2Z
+sequential_4/dense_3/BiasAdd/ReadVariableOp+sequential_4/dense_3/BiasAdd/ReadVariableOp2X
*sequential_4/dense_3/MatMul/ReadVariableOp*sequential_4/dense_3/MatMul/ReadVariableOp2Z
+sequential_4/dense_4/BiasAdd/ReadVariableOp+sequential_4/dense_4/BiasAdd/ReadVariableOp2X
*sequential_4/dense_4/MatMul/ReadVariableOp*sequential_4/dense_4/MatMul/ReadVariableOp2Z
+sequential_4/dense_5/BiasAdd/ReadVariableOp+sequential_4/dense_5/BiasAdd/ReadVariableOp2X
*sequential_4/dense_5/MatMul/ReadVariableOp*sequential_4/dense_5/MatMul/ReadVariableOp:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
`
D__inference_flatten_1_layer_call_and_return_conditional_losses_31528

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"�����  2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:����������2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:����������2

Identity"
identityIdentity:output:0*.
_input_shapes
:��������� :W S
/
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
,__inference_sequential_4_layer_call_fn_31862
sequential_3_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallsequential_3_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_4_layer_call_and_return_conditional_losses_318292
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_3_input
�	
�
B__inference_dense_5_layer_call_and_return_conditional_losses_34065

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
MatMul�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2	
BiasAdda
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:���������2	
Sigmoid�
IdentityIdentitySigmoid:y:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�]
�
G__inference_sequential_4_layer_call_and_return_conditional_losses_33389

inputs+
'conv2d_3_conv2d_readvariableop_resource,
(conv2d_3_biasadd_readvariableop_resource+
'conv2d_4_conv2d_readvariableop_resource,
(conv2d_4_biasadd_readvariableop_resource+
'conv2d_5_conv2d_readvariableop_resource,
(conv2d_5_biasadd_readvariableop_resource*
&dense_3_matmul_readvariableop_resource+
'dense_3_biasadd_readvariableop_resource*
&dense_4_matmul_readvariableop_resource+
'dense_4_biasadd_readvariableop_resource*
&dense_5_matmul_readvariableop_resource+
'dense_5_biasadd_readvariableop_resource
identity��conv2d_3/BiasAdd/ReadVariableOp�conv2d_3/Conv2D/ReadVariableOp�conv2d_4/BiasAdd/ReadVariableOp�conv2d_4/Conv2D/ReadVariableOp�conv2d_5/BiasAdd/ReadVariableOp�conv2d_5/Conv2D/ReadVariableOp�dense_3/BiasAdd/ReadVariableOp�dense_3/MatMul/ReadVariableOp�0dense_3/kernel/Regularizer/Square/ReadVariableOp�dense_4/BiasAdd/ReadVariableOp�dense_4/MatMul/ReadVariableOp�0dense_4/kernel/Regularizer/Square/ReadVariableOp�dense_5/BiasAdd/ReadVariableOp�dense_5/MatMul/ReadVariableOpm
rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2
rescaling_3/Cast/xq
rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rescaling_3/Cast_1/x�
rescaling_3/mulMulinputsrescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/mul�
rescaling_3/addAddV2rescaling_3/mul:z:0rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/add�
conv2d_3/Conv2D/ReadVariableOpReadVariableOp'conv2d_3_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02 
conv2d_3/Conv2D/ReadVariableOp�
conv2d_3/Conv2DConv2Drescaling_3/add:z:0&conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&*
paddingSAME*
strides
2
conv2d_3/Conv2D�
conv2d_3/BiasAdd/ReadVariableOpReadVariableOp(conv2d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv2d_3/BiasAdd/ReadVariableOp�
conv2d_3/BiasAddBiasAddconv2d_3/Conv2D:output:0'conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&2
conv2d_3/BiasAdd{
conv2d_3/ReluReluconv2d_3/BiasAdd:output:0*
T0*/
_output_shapes
:���������&2
conv2d_3/Relu�
dropout_4/IdentityIdentityconv2d_3/Relu:activations:0*
T0*/
_output_shapes
:���������&2
dropout_4/Identity�
max_pooling2d_3/MaxPoolMaxPooldropout_4/Identity:output:0*/
_output_shapes
:���������*
ksize
*
paddingVALID*
strides
2
max_pooling2d_3/MaxPool�
conv2d_4/Conv2D/ReadVariableOpReadVariableOp'conv2d_4_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02 
conv2d_4/Conv2D/ReadVariableOp�
conv2d_4/Conv2DConv2D max_pooling2d_3/MaxPool:output:0&conv2d_4/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������*
paddingSAME*
strides
2
conv2d_4/Conv2D�
conv2d_4/BiasAdd/ReadVariableOpReadVariableOp(conv2d_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv2d_4/BiasAdd/ReadVariableOp�
conv2d_4/BiasAddBiasAddconv2d_4/Conv2D:output:0'conv2d_4/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������2
conv2d_4/BiasAdd{
conv2d_4/ReluReluconv2d_4/BiasAdd:output:0*
T0*/
_output_shapes
:���������2
conv2d_4/Relu�
dropout_5/IdentityIdentityconv2d_4/Relu:activations:0*
T0*/
_output_shapes
:���������2
dropout_5/Identity�
max_pooling2d_4/MaxPoolMaxPooldropout_5/Identity:output:0*/
_output_shapes
:���������	*
ksize
*
paddingVALID*
strides
2
max_pooling2d_4/MaxPool�
conv2d_5/Conv2D/ReadVariableOpReadVariableOp'conv2d_5_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype02 
conv2d_5/Conv2D/ReadVariableOp�
conv2d_5/Conv2DConv2D max_pooling2d_4/MaxPool:output:0&conv2d_5/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 *
paddingSAME*
strides
2
conv2d_5/Conv2D�
conv2d_5/BiasAdd/ReadVariableOpReadVariableOp(conv2d_5_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02!
conv2d_5/BiasAdd/ReadVariableOp�
conv2d_5/BiasAddBiasAddconv2d_5/Conv2D:output:0'conv2d_5/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������	 2
conv2d_5/BiasAdd{
conv2d_5/ReluReluconv2d_5/BiasAdd:output:0*
T0*/
_output_shapes
:���������	 2
conv2d_5/Relu�
max_pooling2d_5/MaxPoolMaxPoolconv2d_5/Relu:activations:0*/
_output_shapes
:��������� *
ksize
*
paddingVALID*
strides
2
max_pooling2d_5/MaxPools
flatten_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"�����  2
flatten_1/Const�
flatten_1/ReshapeReshape max_pooling2d_5/MaxPool:output:0flatten_1/Const:output:0*
T0*(
_output_shapes
:����������2
flatten_1/Reshape�
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype02
dense_3/MatMul/ReadVariableOp�
dense_3/MatMulMatMulflatten_1/Reshape:output:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
dense_3/MatMul�
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype02 
dense_3/BiasAdd/ReadVariableOp�
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������2
dense_3/BiasAddq
dense_3/ReluReludense_3/BiasAdd:output:0*
T0*(
_output_shapes
:����������2
dense_3/Relu�
dropout_6/IdentityIdentitydense_3/Relu:activations:0*
T0*(
_output_shapes
:����������2
dropout_6/Identity�
dense_4/MatMul/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype02
dense_4/MatMul/ReadVariableOp�
dense_4/MatMulMatMuldropout_6/Identity:output:0%dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_4/MatMul�
dense_4/BiasAdd/ReadVariableOpReadVariableOp'dense_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_4/BiasAdd/ReadVariableOp�
dense_4/BiasAddBiasAdddense_4/MatMul:product:0&dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_4/BiasAddp
dense_4/ReluReludense_4/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
dense_4/Relu�
dropout_7/IdentityIdentitydense_4/Relu:activations:0*
T0*'
_output_shapes
:���������2
dropout_7/Identity�
dense_5/MatMul/ReadVariableOpReadVariableOp&dense_5_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
dense_5/MatMul/ReadVariableOp�
dense_5/MatMulMatMuldropout_7/Identity:output:0%dense_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_5/MatMul�
dense_5/BiasAdd/ReadVariableOpReadVariableOp'dense_5_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_5/BiasAdd/ReadVariableOp�
dense_5/BiasAddBiasAdddense_5/MatMul:product:0&dense_5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������2
dense_5/BiasAddy
dense_5/SigmoidSigmoiddense_5/BiasAdd:output:0*
T0*'
_output_shapes
:���������2
dense_5/Sigmoid�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentitydense_5/Sigmoid:y:0 ^conv2d_3/BiasAdd/ReadVariableOp^conv2d_3/Conv2D/ReadVariableOp ^conv2d_4/BiasAdd/ReadVariableOp^conv2d_4/Conv2D/ReadVariableOp ^conv2d_5/BiasAdd/ReadVariableOp^conv2d_5/Conv2D/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp^dense_4/BiasAdd/ReadVariableOp^dense_4/MatMul/ReadVariableOp1^dense_4/kernel/Regularizer/Square/ReadVariableOp^dense_5/BiasAdd/ReadVariableOp^dense_5/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::2B
conv2d_3/BiasAdd/ReadVariableOpconv2d_3/BiasAdd/ReadVariableOp2@
conv2d_3/Conv2D/ReadVariableOpconv2d_3/Conv2D/ReadVariableOp2B
conv2d_4/BiasAdd/ReadVariableOpconv2d_4/BiasAdd/ReadVariableOp2@
conv2d_4/Conv2D/ReadVariableOpconv2d_4/Conv2D/ReadVariableOp2B
conv2d_5/BiasAdd/ReadVariableOpconv2d_5/BiasAdd/ReadVariableOp2@
conv2d_5/Conv2D/ReadVariableOpconv2d_5/Conv2D/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2@
dense_4/BiasAdd/ReadVariableOpdense_4/BiasAdd/ReadVariableOp2>
dense_4/MatMul/ReadVariableOpdense_4/MatMul/ReadVariableOp2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2@
dense_5/BiasAdd/ReadVariableOpdense_5/BiasAdd/ReadVariableOp2>
dense_5/MatMul/ReadVariableOpdense_5/MatMul/ReadVariableOp:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�	
�
,__inference_sequential_5_layer_call_fn_32272
sequential_4_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallsequential_4_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*.
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_5_layer_call_and_return_conditional_losses_322452
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*^
_input_shapesM
K:���������&::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_4_input
�
�
__inference_loss_fn_0_34085=
9dense_3_kernel_regularizer_square_readvariableop_resource
identity��0dense_3/kernel/Regularizer/Square/ReadVariableOp�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_3_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
IdentityIdentity"dense_3/kernel/Regularizer/mul:z:01^dense_3/kernel/Regularizer/Square/ReadVariableOp*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*
_input_shapes
:2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp
�\
�
G__inference_sequential_4_layer_call_and_return_conditional_losses_31702
sequential_3_input
sequential_3_31367
sequential_3_31369
sequential_3_31371
conv2d_3_31400
conv2d_3_31402
conv2d_4_31458
conv2d_4_31460
conv2d_5_31516
conv2d_5_31518
dense_3_31564
dense_3_31566
dense_4_31627
dense_4_31629
dense_5_31684
dense_5_31686
identity�� conv2d_3/StatefulPartitionedCall� conv2d_4/StatefulPartitionedCall� conv2d_5/StatefulPartitionedCall�dense_3/StatefulPartitionedCall�0dense_3/kernel/Regularizer/Square/ReadVariableOp�dense_4/StatefulPartitionedCall�0dense_4/kernel/Regularizer/Square/ReadVariableOp�dense_5/StatefulPartitionedCall�!dropout_4/StatefulPartitionedCall�!dropout_5/StatefulPartitionedCall�!dropout_6/StatefulPartitionedCall�!dropout_7/StatefulPartitionedCall�$sequential_3/StatefulPartitionedCall�
$sequential_3/StatefulPartitionedCallStatefulPartitionedCallsequential_3_inputsequential_3_31367sequential_3_31369sequential_3_31371*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_312932&
$sequential_3/StatefulPartitionedCallm
rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2
rescaling_3/Cast/xq
rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rescaling_3/Cast_1/x�
rescaling_3/mulMul-sequential_3/StatefulPartitionedCall:output:0rescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/mul�
rescaling_3/addAddV2rescaling_3/mul:z:0rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/add�
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCallrescaling_3/add:z:0conv2d_3_31400conv2d_3_31402*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_313892"
 conv2d_3/StatefulPartitionedCall�
!dropout_4/StatefulPartitionedCallStatefulPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_4_layer_call_and_return_conditional_losses_314172#
!dropout_4/StatefulPartitionedCall�
max_pooling2d_3/PartitionedCallPartitionedCall*dropout_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_313172!
max_pooling2d_3/PartitionedCall�
 conv2d_4/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_3/PartitionedCall:output:0conv2d_4_31458conv2d_4_31460*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_314472"
 conv2d_4/StatefulPartitionedCall�
!dropout_5/StatefulPartitionedCallStatefulPartitionedCall)conv2d_4/StatefulPartitionedCall:output:0"^dropout_4/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_5_layer_call_and_return_conditional_losses_314752#
!dropout_5/StatefulPartitionedCall�
max_pooling2d_4/PartitionedCallPartitionedCall*dropout_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_4_layer_call_and_return_conditional_losses_313292!
max_pooling2d_4/PartitionedCall�
 conv2d_5/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_4/PartitionedCall:output:0conv2d_5_31516conv2d_5_31518*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_5_layer_call_and_return_conditional_losses_315052"
 conv2d_5/StatefulPartitionedCall�
max_pooling2d_5/PartitionedCallPartitionedCall)conv2d_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_5_layer_call_and_return_conditional_losses_313412!
max_pooling2d_5/PartitionedCall�
flatten_1/PartitionedCallPartitionedCall(max_pooling2d_5/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_flatten_1_layer_call_and_return_conditional_losses_315282
flatten_1/PartitionedCall�
dense_3/StatefulPartitionedCallStatefulPartitionedCall"flatten_1/PartitionedCall:output:0dense_3_31564dense_3_31566*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_315532!
dense_3/StatefulPartitionedCall�
!dropout_6/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0"^dropout_5/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_6_layer_call_and_return_conditional_losses_315812#
!dropout_6/StatefulPartitionedCall�
dense_4/StatefulPartitionedCallStatefulPartitionedCall*dropout_6/StatefulPartitionedCall:output:0dense_4_31627dense_4_31629*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_4_layer_call_and_return_conditional_losses_316162!
dense_4/StatefulPartitionedCall�
!dropout_7/StatefulPartitionedCallStatefulPartitionedCall(dense_4/StatefulPartitionedCall:output:0"^dropout_6/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_7_layer_call_and_return_conditional_losses_316442#
!dropout_7/StatefulPartitionedCall�
dense_5/StatefulPartitionedCallStatefulPartitionedCall*dropout_7/StatefulPartitionedCall:output:0dense_5_31684dense_5_31686*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_5_layer_call_and_return_conditional_losses_316732!
dense_5/StatefulPartitionedCall�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_31564* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_4_31627*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity(dense_5/StatefulPartitionedCall:output:0!^conv2d_3/StatefulPartitionedCall!^conv2d_4/StatefulPartitionedCall!^conv2d_5/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall1^dense_3/kernel/Regularizer/Square/ReadVariableOp ^dense_4/StatefulPartitionedCall1^dense_4/kernel/Regularizer/Square/ReadVariableOp ^dense_5/StatefulPartitionedCall"^dropout_4/StatefulPartitionedCall"^dropout_5/StatefulPartitionedCall"^dropout_6/StatefulPartitionedCall"^dropout_7/StatefulPartitionedCall%^sequential_3/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2D
 conv2d_4/StatefulPartitionedCall conv2d_4/StatefulPartitionedCall2D
 conv2d_5/StatefulPartitionedCall conv2d_5/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2B
dense_5/StatefulPartitionedCalldense_5/StatefulPartitionedCall2F
!dropout_4/StatefulPartitionedCall!dropout_4/StatefulPartitionedCall2F
!dropout_5/StatefulPartitionedCall!dropout_5/StatefulPartitionedCall2F
!dropout_6/StatefulPartitionedCall!dropout_6/StatefulPartitionedCall2F
!dropout_7/StatefulPartitionedCall!dropout_7/StatefulPartitionedCall2L
$sequential_3/StatefulPartitionedCall$sequential_3/StatefulPartitionedCall:c _
/
_output_shapes
:���������&
,
_user_specified_namesequential_3_input
�

�
C__inference_conv2d_3_layer_call_and_return_conditional_losses_33822

inputs"
conv2d_readvariableop_resource#
biasadd_readvariableop_resource
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&*
paddingSAME*
strides
2
Conv2D�
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp�
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������&2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:���������&2
Relu�
IdentityIdentityRelu:activations:0^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*6
_input_shapes%
#:���������&::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
K
/__inference_max_pooling2d_5_layer_call_fn_31347

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4������������������������������������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_5_layer_call_and_return_conditional_losses_313412
PartitionedCall�
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4������������������������������������2

Identity"
identityIdentity:output:0*I
_input_shapes8
6:4������������������������������������:r n
J
_output_shapes8
6:4������������������������������������
 
_user_specified_nameinputs
��
�
!__inference__traced_restore_34443
file_prefix$
 assignvariableop_conv2d_3_kernel$
 assignvariableop_1_conv2d_3_bias&
"assignvariableop_2_conv2d_4_kernel$
 assignvariableop_3_conv2d_4_bias&
"assignvariableop_4_conv2d_5_kernel$
 assignvariableop_5_conv2d_5_bias%
!assignvariableop_6_dense_3_kernel#
assignvariableop_7_dense_3_bias%
!assignvariableop_8_dense_4_kernel#
assignvariableop_9_dense_4_bias&
"assignvariableop_10_dense_5_kernel$
 assignvariableop_11_dense_5_bias!
assignvariableop_12_adam_iter#
assignvariableop_13_adam_beta_1#
assignvariableop_14_adam_beta_2"
assignvariableop_15_adam_decay*
&assignvariableop_16_adam_learning_rate 
assignvariableop_17_variable"
assignvariableop_18_variable_1"
assignvariableop_19_variable_2"
assignvariableop_20_variable_3
assignvariableop_21_total
assignvariableop_22_count
assignvariableop_23_total_1
assignvariableop_24_count_1.
*assignvariableop_25_adam_conv2d_3_kernel_m,
(assignvariableop_26_adam_conv2d_3_bias_m.
*assignvariableop_27_adam_conv2d_4_kernel_m,
(assignvariableop_28_adam_conv2d_4_bias_m.
*assignvariableop_29_adam_conv2d_5_kernel_m,
(assignvariableop_30_adam_conv2d_5_bias_m-
)assignvariableop_31_adam_dense_3_kernel_m+
'assignvariableop_32_adam_dense_3_bias_m-
)assignvariableop_33_adam_dense_4_kernel_m+
'assignvariableop_34_adam_dense_4_bias_m-
)assignvariableop_35_adam_dense_5_kernel_m+
'assignvariableop_36_adam_dense_5_bias_m.
*assignvariableop_37_adam_conv2d_3_kernel_v,
(assignvariableop_38_adam_conv2d_3_bias_v.
*assignvariableop_39_adam_conv2d_4_kernel_v,
(assignvariableop_40_adam_conv2d_4_bias_v.
*assignvariableop_41_adam_conv2d_5_kernel_v,
(assignvariableop_42_adam_conv2d_5_bias_v-
)assignvariableop_43_adam_dense_3_kernel_v+
'assignvariableop_44_adam_dense_3_bias_v-
)assignvariableop_45_adam_dense_4_kernel_v+
'assignvariableop_46_adam_dense_4_bias_v-
)assignvariableop_47_adam_dense_5_kernel_v+
'assignvariableop_48_adam_dense_5_bias_v
identity_50��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_43�AssignVariableOp_44�AssignVariableOp_45�AssignVariableOp_46�AssignVariableOp_47�AssignVariableOp_48�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9� 
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:2*
dtype0*�
value�B�2B0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUEB>layer_with_weights-0/optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB@layer_with_weights-0/optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB@layer_with_weights-0/optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-0/optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEBGlayer_with_weights-0/optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEBOlayer_with_weights-0/layer-0/layer-0/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEBOlayer_with_weights-0/layer-0/layer-1/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEBOlayer_with_weights-0/layer-0/layer-2/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEBOlayer_with_weights-0/layer-0/layer-3/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEBIlayer_with_weights-0/keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEBIlayer_with_weights-0/keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBIlayer_with_weights-0/keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEBIlayer_with_weights-0/keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/0/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/1/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/2/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/3/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/4/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/5/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/6/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/7/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/8/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/9/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBbtrainable_variables/10/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBbtrainable_variables/11/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/0/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/1/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/2/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/3/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/4/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/5/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/6/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/7/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/8/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBatrainable_variables/9/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBbtrainable_variables/10/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBbtrainable_variables/11/.OPTIMIZER_SLOT/layer_with_weights-0/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_names�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:2*
dtype0*w
valuenBl2B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices�
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�::::::::::::::::::::::::::::::::::::::::::::::::::*@
dtypes6
422					2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

Identity�
AssignVariableOpAssignVariableOp assignvariableop_conv2d_3_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1�
AssignVariableOp_1AssignVariableOp assignvariableop_1_conv2d_3_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2�
AssignVariableOp_2AssignVariableOp"assignvariableop_2_conv2d_4_kernelIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3�
AssignVariableOp_3AssignVariableOp assignvariableop_3_conv2d_4_biasIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4�
AssignVariableOp_4AssignVariableOp"assignvariableop_4_conv2d_5_kernelIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5�
AssignVariableOp_5AssignVariableOp assignvariableop_5_conv2d_5_biasIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6�
AssignVariableOp_6AssignVariableOp!assignvariableop_6_dense_3_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7�
AssignVariableOp_7AssignVariableOpassignvariableop_7_dense_3_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8�
AssignVariableOp_8AssignVariableOp!assignvariableop_8_dense_4_kernelIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9�
AssignVariableOp_9AssignVariableOpassignvariableop_9_dense_4_biasIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10�
AssignVariableOp_10AssignVariableOp"assignvariableop_10_dense_5_kernelIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11�
AssignVariableOp_11AssignVariableOp assignvariableop_11_dense_5_biasIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_12�
AssignVariableOp_12AssignVariableOpassignvariableop_12_adam_iterIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13�
AssignVariableOp_13AssignVariableOpassignvariableop_13_adam_beta_1Identity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14�
AssignVariableOp_14AssignVariableOpassignvariableop_14_adam_beta_2Identity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15�
AssignVariableOp_15AssignVariableOpassignvariableop_15_adam_decayIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16�
AssignVariableOp_16AssignVariableOp&assignvariableop_16_adam_learning_rateIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_17�
AssignVariableOp_17AssignVariableOpassignvariableop_17_variableIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_18�
AssignVariableOp_18AssignVariableOpassignvariableop_18_variable_1Identity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_19�
AssignVariableOp_19AssignVariableOpassignvariableop_19_variable_2Identity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_20�
AssignVariableOp_20AssignVariableOpassignvariableop_20_variable_3Identity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21�
AssignVariableOp_21AssignVariableOpassignvariableop_21_totalIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22�
AssignVariableOp_22AssignVariableOpassignvariableop_22_countIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23�
AssignVariableOp_23AssignVariableOpassignvariableop_23_total_1Identity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24�
AssignVariableOp_24AssignVariableOpassignvariableop_24_count_1Identity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:2
Identity_25�
AssignVariableOp_25AssignVariableOp*assignvariableop_25_adam_conv2d_3_kernel_mIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:2
Identity_26�
AssignVariableOp_26AssignVariableOp(assignvariableop_26_adam_conv2d_3_bias_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:2
Identity_27�
AssignVariableOp_27AssignVariableOp*assignvariableop_27_adam_conv2d_4_kernel_mIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28�
AssignVariableOp_28AssignVariableOp(assignvariableop_28_adam_conv2d_4_bias_mIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29�
AssignVariableOp_29AssignVariableOp*assignvariableop_29_adam_conv2d_5_kernel_mIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30�
AssignVariableOp_30AssignVariableOp(assignvariableop_30_adam_conv2d_5_bias_mIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31�
AssignVariableOp_31AssignVariableOp)assignvariableop_31_adam_dense_3_kernel_mIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32�
AssignVariableOp_32AssignVariableOp'assignvariableop_32_adam_dense_3_bias_mIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33�
AssignVariableOp_33AssignVariableOp)assignvariableop_33_adam_dense_4_kernel_mIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_33n
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:2
Identity_34�
AssignVariableOp_34AssignVariableOp'assignvariableop_34_adam_dense_4_bias_mIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_34n
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:2
Identity_35�
AssignVariableOp_35AssignVariableOp)assignvariableop_35_adam_dense_5_kernel_mIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_35n
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:2
Identity_36�
AssignVariableOp_36AssignVariableOp'assignvariableop_36_adam_dense_5_bias_mIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_36n
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:2
Identity_37�
AssignVariableOp_37AssignVariableOp*assignvariableop_37_adam_conv2d_3_kernel_vIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_37n
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:2
Identity_38�
AssignVariableOp_38AssignVariableOp(assignvariableop_38_adam_conv2d_3_bias_vIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_38n
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:2
Identity_39�
AssignVariableOp_39AssignVariableOp*assignvariableop_39_adam_conv2d_4_kernel_vIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_39n
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:2
Identity_40�
AssignVariableOp_40AssignVariableOp(assignvariableop_40_adam_conv2d_4_bias_vIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_40n
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:2
Identity_41�
AssignVariableOp_41AssignVariableOp*assignvariableop_41_adam_conv2d_5_kernel_vIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_41n
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:2
Identity_42�
AssignVariableOp_42AssignVariableOp(assignvariableop_42_adam_conv2d_5_bias_vIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_42n
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:2
Identity_43�
AssignVariableOp_43AssignVariableOp)assignvariableop_43_adam_dense_3_kernel_vIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_43n
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:2
Identity_44�
AssignVariableOp_44AssignVariableOp'assignvariableop_44_adam_dense_3_bias_vIdentity_44:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_44n
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:2
Identity_45�
AssignVariableOp_45AssignVariableOp)assignvariableop_45_adam_dense_4_kernel_vIdentity_45:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_45n
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:2
Identity_46�
AssignVariableOp_46AssignVariableOp'assignvariableop_46_adam_dense_4_bias_vIdentity_46:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_46n
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:2
Identity_47�
AssignVariableOp_47AssignVariableOp)assignvariableop_47_adam_dense_5_kernel_vIdentity_47:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_47n
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:2
Identity_48�
AssignVariableOp_48AssignVariableOp'assignvariableop_48_adam_dense_5_bias_vIdentity_48:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_489
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp�	
Identity_49Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_49�	
Identity_50IdentityIdentity_49:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*
T0*
_output_shapes
: 2
Identity_50"#
identity_50Identity_50:output:0*�
_input_shapes�
�: :::::::::::::::::::::::::::::::::::::::::::::::::2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�[
�
G__inference_sequential_4_layer_call_and_return_conditional_losses_31829

inputs
sequential_3_31767
sequential_3_31769
sequential_3_31771
conv2d_3_31778
conv2d_3_31780
conv2d_4_31785
conv2d_4_31787
conv2d_5_31792
conv2d_5_31794
dense_3_31799
dense_3_31801
dense_4_31805
dense_4_31807
dense_5_31811
dense_5_31813
identity�� conv2d_3/StatefulPartitionedCall� conv2d_4/StatefulPartitionedCall� conv2d_5/StatefulPartitionedCall�dense_3/StatefulPartitionedCall�0dense_3/kernel/Regularizer/Square/ReadVariableOp�dense_4/StatefulPartitionedCall�0dense_4/kernel/Regularizer/Square/ReadVariableOp�dense_5/StatefulPartitionedCall�!dropout_4/StatefulPartitionedCall�!dropout_5/StatefulPartitionedCall�!dropout_6/StatefulPartitionedCall�!dropout_7/StatefulPartitionedCall�$sequential_3/StatefulPartitionedCall�
$sequential_3/StatefulPartitionedCallStatefulPartitionedCallinputssequential_3_31767sequential_3_31769sequential_3_31771*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_312932&
$sequential_3/StatefulPartitionedCallm
rescaling_3/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *���;2
rescaling_3/Cast/xq
rescaling_3/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rescaling_3/Cast_1/x�
rescaling_3/mulMul-sequential_3/StatefulPartitionedCall:output:0rescaling_3/Cast/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/mul�
rescaling_3/addAddV2rescaling_3/mul:z:0rescaling_3/Cast_1/x:output:0*
T0*/
_output_shapes
:���������&2
rescaling_3/add�
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCallrescaling_3/add:z:0conv2d_3_31778conv2d_3_31780*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_313892"
 conv2d_3/StatefulPartitionedCall�
!dropout_4/StatefulPartitionedCallStatefulPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������&* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_4_layer_call_and_return_conditional_losses_314172#
!dropout_4/StatefulPartitionedCall�
max_pooling2d_3/PartitionedCallPartitionedCall*dropout_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_313172!
max_pooling2d_3/PartitionedCall�
 conv2d_4/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_3/PartitionedCall:output:0conv2d_4_31785conv2d_4_31787*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_314472"
 conv2d_4/StatefulPartitionedCall�
!dropout_5/StatefulPartitionedCallStatefulPartitionedCall)conv2d_4/StatefulPartitionedCall:output:0"^dropout_4/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_5_layer_call_and_return_conditional_losses_314752#
!dropout_5/StatefulPartitionedCall�
max_pooling2d_4/PartitionedCallPartitionedCall*dropout_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_4_layer_call_and_return_conditional_losses_313292!
max_pooling2d_4/PartitionedCall�
 conv2d_5/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_4/PartitionedCall:output:0conv2d_5_31792conv2d_5_31794*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *L
fGRE
C__inference_conv2d_5_layer_call_and_return_conditional_losses_315052"
 conv2d_5/StatefulPartitionedCall�
max_pooling2d_5/PartitionedCallPartitionedCall)conv2d_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_5_layer_call_and_return_conditional_losses_313412!
max_pooling2d_5/PartitionedCall�
flatten_1/PartitionedCallPartitionedCall(max_pooling2d_5/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_flatten_1_layer_call_and_return_conditional_losses_315282
flatten_1/PartitionedCall�
dense_3/StatefulPartitionedCallStatefulPartitionedCall"flatten_1/PartitionedCall:output:0dense_3_31799dense_3_31801*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_315532!
dense_3/StatefulPartitionedCall�
!dropout_6/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0"^dropout_5/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_6_layer_call_and_return_conditional_losses_315812#
!dropout_6/StatefulPartitionedCall�
dense_4/StatefulPartitionedCallStatefulPartitionedCall*dropout_6/StatefulPartitionedCall:output:0dense_4_31805dense_4_31807*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_4_layer_call_and_return_conditional_losses_316162!
dense_4/StatefulPartitionedCall�
!dropout_7/StatefulPartitionedCallStatefulPartitionedCall(dense_4/StatefulPartitionedCall:output:0"^dropout_6/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *M
fHRF
D__inference_dropout_7_layer_call_and_return_conditional_losses_316442#
!dropout_7/StatefulPartitionedCall�
dense_5/StatefulPartitionedCallStatefulPartitionedCall*dropout_7/StatefulPartitionedCall:output:0dense_5_31811dense_5_31813*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_5_layer_call_and_return_conditional_losses_316732!
dense_5/StatefulPartitionedCall�
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_31799* 
_output_shapes
:
��*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp�
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
��2#
!dense_3/kernel/Regularizer/Square�
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const�
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum�
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_3/kernel/Regularizer/mul/x�
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul�
0dense_4/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_4_31805*
_output_shapes
:	�*
dtype022
0dense_4/kernel/Regularizer/Square/ReadVariableOp�
!dense_4/kernel/Regularizer/SquareSquare8dense_4/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	�2#
!dense_4/kernel/Regularizer/Square�
 dense_4/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_4/kernel/Regularizer/Const�
dense_4/kernel/Regularizer/SumSum%dense_4/kernel/Regularizer/Square:y:0)dense_4/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/Sum�
 dense_4/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o�:2"
 dense_4/kernel/Regularizer/mul/x�
dense_4/kernel/Regularizer/mulMul)dense_4/kernel/Regularizer/mul/x:output:0'dense_4/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_4/kernel/Regularizer/mul�
IdentityIdentity(dense_5/StatefulPartitionedCall:output:0!^conv2d_3/StatefulPartitionedCall!^conv2d_4/StatefulPartitionedCall!^conv2d_5/StatefulPartitionedCall ^dense_3/StatefulPartitionedCall1^dense_3/kernel/Regularizer/Square/ReadVariableOp ^dense_4/StatefulPartitionedCall1^dense_4/kernel/Regularizer/Square/ReadVariableOp ^dense_5/StatefulPartitionedCall"^dropout_4/StatefulPartitionedCall"^dropout_5/StatefulPartitionedCall"^dropout_6/StatefulPartitionedCall"^dropout_7/StatefulPartitionedCall%^sequential_3/StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*j
_input_shapesY
W:���������&:::::::::::::::2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2D
 conv2d_4/StatefulPartitionedCall conv2d_4/StatefulPartitionedCall2D
 conv2d_5/StatefulPartitionedCall conv2d_5/StatefulPartitionedCall2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall2d
0dense_4/kernel/Regularizer/Square/ReadVariableOp0dense_4/kernel/Regularizer/Square/ReadVariableOp2B
dense_5/StatefulPartitionedCalldense_5/StatefulPartitionedCall2F
!dropout_4/StatefulPartitionedCall!dropout_4/StatefulPartitionedCall2F
!dropout_5/StatefulPartitionedCall!dropout_5/StatefulPartitionedCall2F
!dropout_6/StatefulPartitionedCall!dropout_6/StatefulPartitionedCall2F
!dropout_7/StatefulPartitionedCall!dropout_7/StatefulPartitionedCall2L
$sequential_3/StatefulPartitionedCall$sequential_3/StatefulPartitionedCall:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
K
/__inference_max_pooling2d_3_layer_call_fn_31323

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4������������������������������������* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_313172
PartitionedCall�
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4������������������������������������2

Identity"
identityIdentity:output:0*I
_input_shapes8
6:4������������������������������������:r n
J
_output_shapes8
6:4������������������������������������
 
_user_specified_nameinputs
�
c
D__inference_dropout_4_layer_call_and_return_conditional_losses_31417

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *�8�?2
dropout/Const{
dropout/MulMulinputsdropout/Const:output:0*
T0*/
_output_shapes
:���������&2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*/
_output_shapes
:���������&*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���=2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*/
_output_shapes
:���������&2
dropout/GreaterEqual�
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*/
_output_shapes
:���������&2
dropout/Cast�
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*/
_output_shapes
:���������&2
dropout/Mul_1m
IdentityIdentitydropout/Mul_1:z:0*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&:W S
/
_output_shapes
:���������&
 
_user_specified_nameinputs
�
p
G__inference_sequential_3_layer_call_and_return_conditional_losses_30962
random_flip_1_input
identityo
IdentityIdentityrandom_flip_1_input*
T0*/
_output_shapes
:���������&2

Identity"
identityIdentity:output:0*.
_input_shapes
:���������&:d `
/
_output_shapes
:���������&
-
_user_specified_namerandom_flip_1_input
�
c
D__inference_dropout_7_layer_call_and_return_conditional_losses_34039

inputs
identity�c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:���������2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:���������*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *���>2
dropout/GreaterEqual/y�
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:���������2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:���������2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:���������2
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*&
_input_shapes
:���������:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
|
'__inference_dense_4_layer_call_fn_34027

inputs
unknown
	unknown_0
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_4_layer_call_and_return_conditional_losses_316162
StatefulPartitionedCall�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������2

Identity"
identityIdentity:output:0*/
_input_shapes
:����������::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs"�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
Y
sequential_4_inputC
$serving_default_sequential_4_input:0���������&=
	softmax_10
StatefulPartitionedCall:0���������tensorflow/serving/predict:�
�t
layer_with_weights-0
layer-0
layer-1
regularization_losses
trainable_variables
	variables
	keras_api

signatures
+�&call_and_return_all_conditional_losses
�__call__
�_default_save_signature"�r
_tf_keras_sequential�r{"class_name": "Sequential", "name": "sequential_5", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "sequential_5", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "sequential_4_input"}}, {"class_name": "Sequential", "config": {"name": "sequential_4", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "sequential_3_input"}}, {"class_name": "Sequential", "config": {"name": "sequential_3", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "random_flip_1_input"}}, {"class_name": "RandomFlip", "config": {"name": "random_flip_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "mode": "horizontal", "seed": null}}, {"class_name": "RandomRotation", "config": {"name": "random_rotation_1", "trainable": true, "dtype": "float32", "factor": 0.0, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomTranslation", "config": {"name": "random_translation_1", "trainable": true, "dtype": "float32", "height_factor": {"class_name": "__tuple__", "items": [0, 0.15]}, "width_factor": 0, "fill_mode": "nearest", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomZoom", "config": {"name": "random_zoom_1", "trainable": true, "dtype": "float32", "height_factor": 0.2, "width_factor": null, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}]}}, {"class_name": "Rescaling", "config": {"name": "rescaling_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "scale": 0.00392156862745098, "offset": 0.0}}, {"class_name": "Conv2D", "config": {"name": "conv2d_3", "trainable": true, "dtype": "float32", "filters": 8, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_4", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_3", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Conv2D", "config": {"name": "conv2d_4", "trainable": true, "dtype": "float32", "filters": 16, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_5", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_4", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Conv2D", "config": {"name": "conv2d_5", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_5", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Flatten", "config": {"name": "flatten_1", "trainable": true, "dtype": "float32", "data_format": "channels_last"}}, {"class_name": "Dense", "config": {"name": "dense_3", "trainable": true, "dtype": "float32", "units": 250, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_6", "trainable": true, "dtype": "float32", "rate": 0.5, "noise_shape": null, "seed": null}}, {"class_name": "Dense", "config": {"name": "dense_4", "trainable": true, "dtype": "float32", "units": 20, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_7", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}, {"class_name": "Dense", "config": {"name": "dense_5", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}}, {"class_name": "Softmax", "config": {"name": "softmax_1", "trainable": true, "dtype": "float32", "axis": -1}}]}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 4, "max_ndim": null, "min_ndim": null, "axes": {}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 38, 28, 1]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Sequential", "config": {"name": "sequential_5", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "sequential_4_input"}}, {"class_name": "Sequential", "config": {"name": "sequential_4", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "sequential_3_input"}}, {"class_name": "Sequential", "config": {"name": "sequential_3", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "random_flip_1_input"}}, {"class_name": "RandomFlip", "config": {"name": "random_flip_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "mode": "horizontal", "seed": null}}, {"class_name": "RandomRotation", "config": {"name": "random_rotation_1", "trainable": true, "dtype": "float32", "factor": 0.0, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomTranslation", "config": {"name": "random_translation_1", "trainable": true, "dtype": "float32", "height_factor": {"class_name": "__tuple__", "items": [0, 0.15]}, "width_factor": 0, "fill_mode": "nearest", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomZoom", "config": {"name": "random_zoom_1", "trainable": true, "dtype": "float32", "height_factor": 0.2, "width_factor": null, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}]}}, {"class_name": "Rescaling", "config": {"name": "rescaling_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "scale": 0.00392156862745098, "offset": 0.0}}, {"class_name": "Conv2D", "config": {"name": "conv2d_3", "trainable": true, "dtype": "float32", "filters": 8, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_4", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_3", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Conv2D", "config": {"name": "conv2d_4", "trainable": true, "dtype": "float32", "filters": 16, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_5", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_4", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Conv2D", "config": {"name": "conv2d_5", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_5", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Flatten", "config": {"name": "flatten_1", "trainable": true, "dtype": "float32", "data_format": "channels_last"}}, {"class_name": "Dense", "config": {"name": "dense_3", "trainable": true, "dtype": "float32", "units": 250, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_6", "trainable": true, "dtype": "float32", "rate": 0.5, "noise_shape": null, "seed": null}}, {"class_name": "Dense", "config": {"name": "dense_4", "trainable": true, "dtype": "float32", "units": 20, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_7", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}, {"class_name": "Dense", "config": {"name": "dense_5", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}}, {"class_name": "Softmax", "config": {"name": "softmax_1", "trainable": true, "dtype": "float32", "axis": -1}}]}}}
�t
layer-0
	layer-1

layer_with_weights-0

layer-2
layer-3
layer-4
layer_with_weights-1
layer-5
layer-6
layer-7
layer_with_weights-2
layer-8
layer-9
layer-10
layer_with_weights-3
layer-11
layer-12
layer_with_weights-4
layer-13
layer-14
layer_with_weights-5
layer-15
	optimizer
regularization_losses
trainable_variables
	variables
	keras_api
+�&call_and_return_all_conditional_losses
�__call__"�p
_tf_keras_sequential�p{"class_name": "Sequential", "name": "sequential_4", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "sequential_4", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "sequential_3_input"}}, {"class_name": "Sequential", "config": {"name": "sequential_3", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "random_flip_1_input"}}, {"class_name": "RandomFlip", "config": {"name": "random_flip_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "mode": "horizontal", "seed": null}}, {"class_name": "RandomRotation", "config": {"name": "random_rotation_1", "trainable": true, "dtype": "float32", "factor": 0.0, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomTranslation", "config": {"name": "random_translation_1", "trainable": true, "dtype": "float32", "height_factor": {"class_name": "__tuple__", "items": [0, 0.15]}, "width_factor": 0, "fill_mode": "nearest", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomZoom", "config": {"name": "random_zoom_1", "trainable": true, "dtype": "float32", "height_factor": 0.2, "width_factor": null, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}]}}, {"class_name": "Rescaling", "config": {"name": "rescaling_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "scale": 0.00392156862745098, "offset": 0.0}}, {"class_name": "Conv2D", "config": {"name": "conv2d_3", "trainable": true, "dtype": "float32", "filters": 8, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_4", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_3", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Conv2D", "config": {"name": "conv2d_4", "trainable": true, "dtype": "float32", "filters": 16, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_5", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_4", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Conv2D", "config": {"name": "conv2d_5", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_5", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Flatten", "config": {"name": "flatten_1", "trainable": true, "dtype": "float32", "data_format": "channels_last"}}, {"class_name": "Dense", "config": {"name": "dense_3", "trainable": true, "dtype": "float32", "units": 250, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_6", "trainable": true, "dtype": "float32", "rate": 0.5, "noise_shape": null, "seed": null}}, {"class_name": "Dense", "config": {"name": "dense_4", "trainable": true, "dtype": "float32", "units": 20, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_7", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}, {"class_name": "Dense", "config": {"name": "dense_5", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 4, "max_ndim": null, "min_ndim": null, "axes": {}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 38, 28, 1]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Sequential", "config": {"name": "sequential_4", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "sequential_3_input"}}, {"class_name": "Sequential", "config": {"name": "sequential_3", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "random_flip_1_input"}}, {"class_name": "RandomFlip", "config": {"name": "random_flip_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "mode": "horizontal", "seed": null}}, {"class_name": "RandomRotation", "config": {"name": "random_rotation_1", "trainable": true, "dtype": "float32", "factor": 0.0, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomTranslation", "config": {"name": "random_translation_1", "trainable": true, "dtype": "float32", "height_factor": {"class_name": "__tuple__", "items": [0, 0.15]}, "width_factor": 0, "fill_mode": "nearest", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomZoom", "config": {"name": "random_zoom_1", "trainable": true, "dtype": "float32", "height_factor": 0.2, "width_factor": null, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}]}}, {"class_name": "Rescaling", "config": {"name": "rescaling_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "scale": 0.00392156862745098, "offset": 0.0}}, {"class_name": "Conv2D", "config": {"name": "conv2d_3", "trainable": true, "dtype": "float32", "filters": 8, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_4", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_3", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Conv2D", "config": {"name": "conv2d_4", "trainable": true, "dtype": "float32", "filters": 16, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_5", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_4", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Conv2D", "config": {"name": "conv2d_5", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "MaxPooling2D", "config": {"name": "max_pooling2d_5", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}}, {"class_name": "Flatten", "config": {"name": "flatten_1", "trainable": true, "dtype": "float32", "data_format": "channels_last"}}, {"class_name": "Dense", "config": {"name": "dense_3", "trainable": true, "dtype": "float32", "units": 250, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_6", "trainable": true, "dtype": "float32", "rate": 0.5, "noise_shape": null, "seed": null}}, {"class_name": "Dense", "config": {"name": "dense_4", "trainable": true, "dtype": "float32", "units": 20, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dropout", "config": {"name": "dropout_7", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}, {"class_name": "Dense", "config": {"name": "dense_5", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}}, "training_config": {"loss": "binary_crossentropy", "metrics": [[{"class_name": "MeanMetricWrapper", "config": {"name": "accuracy", "dtype": "float32", "fn": "binary_accuracy"}}]], "weighted_metrics": null, "loss_weights": null, "optimizer_config": {"class_name": "Adam", "config": {"name": "Adam", "learning_rate": 0.0010000000474974513, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07, "amsgrad": false}}}}
�
regularization_losses
trainable_variables
	variables
 	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Softmax", "name": "softmax_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "softmax_1", "trainable": true, "dtype": "float32", "axis": -1}}
 "
trackable_list_wrapper
v
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11"
trackable_list_wrapper
v
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11"
trackable_list_wrapper
�
regularization_losses
-layer_metrics
.layer_regularization_losses
/metrics
trainable_variables

0layers
1non_trainable_variables
	variables
�__call__
�_default_save_signature
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
-
�serving_default"
signature_map
�
2layer-0
3layer-1
4layer-2
5layer-3
6regularization_losses
7trainable_variables
8	variables
9	keras_api
+�&call_and_return_all_conditional_losses
�__call__"�
_tf_keras_sequential�{"class_name": "Sequential", "name": "sequential_3", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "sequential_3", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "random_flip_1_input"}}, {"class_name": "RandomFlip", "config": {"name": "random_flip_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "mode": "horizontal", "seed": null}}, {"class_name": "RandomRotation", "config": {"name": "random_rotation_1", "trainable": true, "dtype": "float32", "factor": 0.0, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomTranslation", "config": {"name": "random_translation_1", "trainable": true, "dtype": "float32", "height_factor": {"class_name": "__tuple__", "items": [0, 0.15]}, "width_factor": 0, "fill_mode": "nearest", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomZoom", "config": {"name": "random_zoom_1", "trainable": true, "dtype": "float32", "height_factor": 0.2, "width_factor": null, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}]}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 4, "max_ndim": null, "min_ndim": null, "axes": {}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 38, 28, 1]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Sequential", "config": {"name": "sequential_3", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "random_flip_1_input"}}, {"class_name": "RandomFlip", "config": {"name": "random_flip_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "mode": "horizontal", "seed": null}}, {"class_name": "RandomRotation", "config": {"name": "random_rotation_1", "trainable": true, "dtype": "float32", "factor": 0.0, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomTranslation", "config": {"name": "random_translation_1", "trainable": true, "dtype": "float32", "height_factor": {"class_name": "__tuple__", "items": [0, 0.15]}, "width_factor": 0, "fill_mode": "nearest", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}, {"class_name": "RandomZoom", "config": {"name": "random_zoom_1", "trainable": true, "dtype": "float32", "height_factor": 0.2, "width_factor": null, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}]}}}
�
:	keras_api"�
_tf_keras_layer�{"class_name": "Rescaling", "name": "rescaling_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "stateful": false, "must_restore_from_config": true, "config": {"name": "rescaling_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "scale": 0.00392156862745098, "offset": 0.0}}
�	

!kernel
"bias
;regularization_losses
<trainable_variables
=	variables
>	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Conv2D", "name": "conv2d_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv2d_3", "trainable": true, "dtype": "float32", "filters": 8, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 4, "axes": {"-1": 1}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 38, 28, 1]}}
�
?regularization_losses
@trainable_variables
A	variables
B	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dropout", "name": "dropout_4", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dropout_4", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}
�
Cregularization_losses
Dtrainable_variables
E	variables
F	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "MaxPooling2D", "name": "max_pooling2d_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "max_pooling2d_3", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 4, "max_ndim": null, "min_ndim": null, "axes": {}}}}
�	

#kernel
$bias
Gregularization_losses
Htrainable_variables
I	variables
J	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Conv2D", "name": "conv2d_4", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv2d_4", "trainable": true, "dtype": "float32", "filters": 16, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 4, "axes": {"-1": 8}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 19, 14, 8]}}
�
Kregularization_losses
Ltrainable_variables
M	variables
N	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dropout", "name": "dropout_5", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dropout_5", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}
�
Oregularization_losses
Ptrainable_variables
Q	variables
R	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "MaxPooling2D", "name": "max_pooling2d_4", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "max_pooling2d_4", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 4, "max_ndim": null, "min_ndim": null, "axes": {}}}}
�	

%kernel
&bias
Sregularization_losses
Ttrainable_variables
U	variables
V	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Conv2D", "name": "conv2d_5", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv2d_5", "trainable": true, "dtype": "float32", "filters": 32, "kernel_size": {"class_name": "__tuple__", "items": [3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1]}, "padding": "same", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1]}, "groups": 1, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 4, "axes": {"-1": 16}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 9, 7, 16]}}
�
Wregularization_losses
Xtrainable_variables
Y	variables
Z	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "MaxPooling2D", "name": "max_pooling2d_5", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "max_pooling2d_5", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [2, 2]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [2, 2]}, "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 4, "max_ndim": null, "min_ndim": null, "axes": {}}}}
�
[regularization_losses
\trainable_variables
]	variables
^	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Flatten", "name": "flatten_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "flatten_1", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}}}
�

'kernel
(bias
_regularization_losses
`trainable_variables
a	variables
b	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_3", "trainable": true, "dtype": "float32", "units": 250, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 384}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 384]}}
�
cregularization_losses
dtrainable_variables
e	variables
f	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dropout", "name": "dropout_6", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dropout_6", "trainable": true, "dtype": "float32", "rate": 0.5, "noise_shape": null, "seed": null}}
�

)kernel
*bias
gregularization_losses
htrainable_variables
i	variables
j	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_4", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_4", "trainable": true, "dtype": "float32", "units": 20, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 250}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 250]}}
�
kregularization_losses
ltrainable_variables
m	variables
n	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dropout", "name": "dropout_7", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dropout_7", "trainable": true, "dtype": "float32", "rate": 0.3, "noise_shape": null, "seed": null}}
�

+kernel
,bias
oregularization_losses
ptrainable_variables
q	variables
r	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "dense_5", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_5", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 20}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 20]}}
�
siter

tbeta_1

ubeta_2
	vdecay
wlearning_rate!m�"m�#m�$m�%m�&m�'m�(m�)m�*m�+m�,m�!v�"v�#v�$v�%v�&v�'v�(v�)v�*v�+v�,v�"
	optimizer
0
�0
�1"
trackable_list_wrapper
v
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11"
trackable_list_wrapper
v
!0
"1
#2
$3
%4
&5
'6
(7
)8
*9
+10
,11"
trackable_list_wrapper
�
regularization_losses
xlayer_metrics
ylayer_regularization_losses
zmetrics
trainable_variables

{layers
|non_trainable_variables
	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
regularization_losses
}layer_metrics
~layer_regularization_losses
metrics
trainable_variables
�layers
�non_trainable_variables
	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
):'2conv2d_3/kernel
:2conv2d_3/bias
):'2conv2d_4/kernel
:2conv2d_4/bias
):' 2conv2d_5/kernel
: 2conv2d_5/bias
": 
��2dense_3/kernel
:�2dense_3/bias
!:	�2dense_4/kernel
:2dense_4/bias
 :2dense_5/kernel
:2dense_5/bias
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
	�_rng
�	keras_api"�
_tf_keras_layer�{"class_name": "RandomFlip", "name": "random_flip_1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "stateful": false, "must_restore_from_config": true, "config": {"name": "random_flip_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 38, 28, 1]}, "dtype": "float32", "mode": "horizontal", "seed": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 4, "max_ndim": null, "min_ndim": null, "axes": {}}}}
�
	�_rng
�	keras_api"�
_tf_keras_layer�{"class_name": "RandomRotation", "name": "random_rotation_1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": true, "config": {"name": "random_rotation_1", "trainable": true, "dtype": "float32", "factor": 0.0, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}
�
	�_rng
�	keras_api"�
_tf_keras_layer�{"class_name": "RandomTranslation", "name": "random_translation_1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": true, "config": {"name": "random_translation_1", "trainable": true, "dtype": "float32", "height_factor": {"class_name": "__tuple__", "items": [0, 0.15]}, "width_factor": 0, "fill_mode": "nearest", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}
�
	�_rng
�	keras_api"�
_tf_keras_layer�{"class_name": "RandomZoom", "name": "random_zoom_1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": true, "config": {"name": "random_zoom_1", "trainable": true, "dtype": "float32", "height_factor": 0.2, "width_factor": null, "fill_mode": "reflect", "fill_value": 0.0, "interpolation": "bilinear", "seed": null}}
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
6regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
7trainable_variables
�layers
�non_trainable_variables
8	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
"
_generic_user_object
 "
trackable_list_wrapper
.
!0
"1"
trackable_list_wrapper
.
!0
"1"
trackable_list_wrapper
�
;regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
<trainable_variables
�layers
�non_trainable_variables
=	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
?regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
@trainable_variables
�layers
�non_trainable_variables
A	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
Cregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Dtrainable_variables
�layers
�non_trainable_variables
E	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
.
#0
$1"
trackable_list_wrapper
�
Gregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Htrainable_variables
�layers
�non_trainable_variables
I	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
Kregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Ltrainable_variables
�layers
�non_trainable_variables
M	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
Oregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Ptrainable_variables
�layers
�non_trainable_variables
Q	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
�
Sregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Ttrainable_variables
�layers
�non_trainable_variables
U	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
Wregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
Xtrainable_variables
�layers
�non_trainable_variables
Y	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
[regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
\trainable_variables
�layers
�non_trainable_variables
]	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
(
�0"
trackable_list_wrapper
.
'0
(1"
trackable_list_wrapper
.
'0
(1"
trackable_list_wrapper
�
_regularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
`trainable_variables
�layers
�non_trainable_variables
a	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
cregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
dtrainable_variables
�layers
�non_trainable_variables
e	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
(
�0"
trackable_list_wrapper
.
)0
*1"
trackable_list_wrapper
.
)0
*1"
trackable_list_wrapper
�
gregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
htrainable_variables
�layers
�non_trainable_variables
i	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
kregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
ltrainable_variables
�layers
�non_trainable_variables
m	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
�
oregularization_losses
�layer_metrics
 �layer_regularization_losses
�metrics
ptrainable_variables
�layers
�non_trainable_variables
q	variables
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
�
0
	1

2
3
4
5
6
7
8
9
10
11
12
13
14
15"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
/
�
_state_var"
_generic_user_object
"
_generic_user_object
/
�
_state_var"
_generic_user_object
"
_generic_user_object
/
�
_state_var"
_generic_user_object
"
_generic_user_object
/
�
_state_var"
_generic_user_object
"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
<
20
31
42
53"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
(
�0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
(
�0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�

�total

�count
�	variables
�	keras_api"�
_tf_keras_metricj{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}}
�

�total

�count
�
_fn_kwargs
�	variables
�	keras_api"�
_tf_keras_metric�{"class_name": "MeanMetricWrapper", "name": "accuracy", "dtype": "float32", "config": {"name": "accuracy", "dtype": "float32", "fn": "binary_accuracy"}}
:	2Variable
:	2Variable
:	2Variable
:	2Variable
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
.:,2Adam/conv2d_3/kernel/m
 :2Adam/conv2d_3/bias/m
.:,2Adam/conv2d_4/kernel/m
 :2Adam/conv2d_4/bias/m
.:, 2Adam/conv2d_5/kernel/m
 : 2Adam/conv2d_5/bias/m
':%
��2Adam/dense_3/kernel/m
 :�2Adam/dense_3/bias/m
&:$	�2Adam/dense_4/kernel/m
:2Adam/dense_4/bias/m
%:#2Adam/dense_5/kernel/m
:2Adam/dense_5/bias/m
.:,2Adam/conv2d_3/kernel/v
 :2Adam/conv2d_3/bias/v
.:,2Adam/conv2d_4/kernel/v
 :2Adam/conv2d_4/bias/v
.:, 2Adam/conv2d_5/kernel/v
 : 2Adam/conv2d_5/bias/v
':%
��2Adam/dense_3/kernel/v
 :�2Adam/dense_3/bias/v
&:$	�2Adam/dense_4/kernel/v
:2Adam/dense_4/bias/v
%:#2Adam/dense_5/kernel/v
:2Adam/dense_5/bias/v
�2�
G__inference_sequential_5_layer_call_and_return_conditional_losses_32117
G__inference_sequential_5_layer_call_and_return_conditional_losses_32739
G__inference_sequential_5_layer_call_and_return_conditional_losses_32075
G__inference_sequential_5_layer_call_and_return_conditional_losses_32811�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
,__inference_sequential_5_layer_call_fn_32875
,__inference_sequential_5_layer_call_fn_32272
,__inference_sequential_5_layer_call_fn_32201
,__inference_sequential_5_layer_call_fn_32846�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
 __inference__wrapped_model_30629�
���
FullArgSpec
args� 
varargsjargs
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *9�6
4�1
sequential_4_input���������&
�2�
G__inference_sequential_4_layer_call_and_return_conditional_losses_33318
G__inference_sequential_4_layer_call_and_return_conditional_losses_31761
G__inference_sequential_4_layer_call_and_return_conditional_losses_33389
G__inference_sequential_4_layer_call_and_return_conditional_losses_31702�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
,__inference_sequential_4_layer_call_fn_33453
,__inference_sequential_4_layer_call_fn_31950
,__inference_sequential_4_layer_call_fn_31862
,__inference_sequential_4_layer_call_fn_33424�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
)__inference_softmax_1_layer_call_fn_33463�
���
FullArgSpec%
args�
jself
jinputs
jmask
varargs
 
varkw
 
defaults�

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_softmax_1_layer_call_and_return_conditional_losses_33458�
���
FullArgSpec%
args�
jself
jinputs
jmask
varargs
 
varkw
 
defaults�

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
#__inference_signature_wrapper_32315sequential_4_input"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
G__inference_sequential_3_layer_call_and_return_conditional_losses_33795
G__inference_sequential_3_layer_call_and_return_conditional_losses_33791
G__inference_sequential_3_layer_call_and_return_conditional_losses_30958
G__inference_sequential_3_layer_call_and_return_conditional_losses_30962�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
,__inference_sequential_3_layer_call_fn_33806
,__inference_sequential_3_layer_call_fn_31311
,__inference_sequential_3_layer_call_fn_31302
,__inference_sequential_3_layer_call_fn_33811�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
(__inference_conv2d_3_layer_call_fn_33831�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_conv2d_3_layer_call_and_return_conditional_losses_33822�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dropout_4_layer_call_fn_33858
)__inference_dropout_4_layer_call_fn_33853�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
D__inference_dropout_4_layer_call_and_return_conditional_losses_33843
D__inference_dropout_4_layer_call_and_return_conditional_losses_33848�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
/__inference_max_pooling2d_3_layer_call_fn_31323�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *@�=
;�84������������������������������������
�2�
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_31317�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *@�=
;�84������������������������������������
�2�
(__inference_conv2d_4_layer_call_fn_33878�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_conv2d_4_layer_call_and_return_conditional_losses_33869�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dropout_5_layer_call_fn_33905
)__inference_dropout_5_layer_call_fn_33900�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
D__inference_dropout_5_layer_call_and_return_conditional_losses_33890
D__inference_dropout_5_layer_call_and_return_conditional_losses_33895�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
/__inference_max_pooling2d_4_layer_call_fn_31335�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *@�=
;�84������������������������������������
�2�
J__inference_max_pooling2d_4_layer_call_and_return_conditional_losses_31329�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *@�=
;�84������������������������������������
�2�
(__inference_conv2d_5_layer_call_fn_33925�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
C__inference_conv2d_5_layer_call_and_return_conditional_losses_33916�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
/__inference_max_pooling2d_5_layer_call_fn_31347�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *@�=
;�84������������������������������������
�2�
J__inference_max_pooling2d_5_layer_call_and_return_conditional_losses_31341�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *@�=
;�84������������������������������������
�2�
)__inference_flatten_1_layer_call_fn_33936�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
D__inference_flatten_1_layer_call_and_return_conditional_losses_33931�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
'__inference_dense_3_layer_call_fn_33968�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
B__inference_dense_3_layer_call_and_return_conditional_losses_33959�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dropout_6_layer_call_fn_33990
)__inference_dropout_6_layer_call_fn_33995�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
D__inference_dropout_6_layer_call_and_return_conditional_losses_33980
D__inference_dropout_6_layer_call_and_return_conditional_losses_33985�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
'__inference_dense_4_layer_call_fn_34027�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
B__inference_dense_4_layer_call_and_return_conditional_losses_34018�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
)__inference_dropout_7_layer_call_fn_34054
)__inference_dropout_7_layer_call_fn_34049�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
D__inference_dropout_7_layer_call_and_return_conditional_losses_34039
D__inference_dropout_7_layer_call_and_return_conditional_losses_34044�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
'__inference_dense_5_layer_call_fn_34074�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
B__inference_dense_5_layer_call_and_return_conditional_losses_34065�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
__inference_loss_fn_0_34085�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�2�
__inference_loss_fn_1_34096�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� �
 __inference__wrapped_model_30629�!"#$%&'()*+,C�@
9�6
4�1
sequential_4_input���������&
� "5�2
0
	softmax_1#� 
	softmax_1����������
C__inference_conv2d_3_layer_call_and_return_conditional_losses_33822l!"7�4
-�*
(�%
inputs���������&
� "-�*
#� 
0���������&
� �
(__inference_conv2d_3_layer_call_fn_33831_!"7�4
-�*
(�%
inputs���������&
� " ����������&�
C__inference_conv2d_4_layer_call_and_return_conditional_losses_33869l#$7�4
-�*
(�%
inputs���������
� "-�*
#� 
0���������
� �
(__inference_conv2d_4_layer_call_fn_33878_#$7�4
-�*
(�%
inputs���������
� " �����������
C__inference_conv2d_5_layer_call_and_return_conditional_losses_33916l%&7�4
-�*
(�%
inputs���������	
� "-�*
#� 
0���������	 
� �
(__inference_conv2d_5_layer_call_fn_33925_%&7�4
-�*
(�%
inputs���������	
� " ����������	 �
B__inference_dense_3_layer_call_and_return_conditional_losses_33959^'(0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� |
'__inference_dense_3_layer_call_fn_33968Q'(0�-
&�#
!�
inputs����������
� "������������
B__inference_dense_4_layer_call_and_return_conditional_losses_34018])*0�-
&�#
!�
inputs����������
� "%�"
�
0���������
� {
'__inference_dense_4_layer_call_fn_34027P)*0�-
&�#
!�
inputs����������
� "�����������
B__inference_dense_5_layer_call_and_return_conditional_losses_34065\+,/�,
%�"
 �
inputs���������
� "%�"
�
0���������
� z
'__inference_dense_5_layer_call_fn_34074O+,/�,
%�"
 �
inputs���������
� "�����������
D__inference_dropout_4_layer_call_and_return_conditional_losses_33843l;�8
1�.
(�%
inputs���������&
p
� "-�*
#� 
0���������&
� �
D__inference_dropout_4_layer_call_and_return_conditional_losses_33848l;�8
1�.
(�%
inputs���������&
p 
� "-�*
#� 
0���������&
� �
)__inference_dropout_4_layer_call_fn_33853_;�8
1�.
(�%
inputs���������&
p
� " ����������&�
)__inference_dropout_4_layer_call_fn_33858_;�8
1�.
(�%
inputs���������&
p 
� " ����������&�
D__inference_dropout_5_layer_call_and_return_conditional_losses_33890l;�8
1�.
(�%
inputs���������
p
� "-�*
#� 
0���������
� �
D__inference_dropout_5_layer_call_and_return_conditional_losses_33895l;�8
1�.
(�%
inputs���������
p 
� "-�*
#� 
0���������
� �
)__inference_dropout_5_layer_call_fn_33900_;�8
1�.
(�%
inputs���������
p
� " �����������
)__inference_dropout_5_layer_call_fn_33905_;�8
1�.
(�%
inputs���������
p 
� " �����������
D__inference_dropout_6_layer_call_and_return_conditional_losses_33980^4�1
*�'
!�
inputs����������
p
� "&�#
�
0����������
� �
D__inference_dropout_6_layer_call_and_return_conditional_losses_33985^4�1
*�'
!�
inputs����������
p 
� "&�#
�
0����������
� ~
)__inference_dropout_6_layer_call_fn_33990Q4�1
*�'
!�
inputs����������
p
� "�����������~
)__inference_dropout_6_layer_call_fn_33995Q4�1
*�'
!�
inputs����������
p 
� "������������
D__inference_dropout_7_layer_call_and_return_conditional_losses_34039\3�0
)�&
 �
inputs���������
p
� "%�"
�
0���������
� �
D__inference_dropout_7_layer_call_and_return_conditional_losses_34044\3�0
)�&
 �
inputs���������
p 
� "%�"
�
0���������
� |
)__inference_dropout_7_layer_call_fn_34049O3�0
)�&
 �
inputs���������
p
� "����������|
)__inference_dropout_7_layer_call_fn_34054O3�0
)�&
 �
inputs���������
p 
� "�����������
D__inference_flatten_1_layer_call_and_return_conditional_losses_33931a7�4
-�*
(�%
inputs��������� 
� "&�#
�
0����������
� �
)__inference_flatten_1_layer_call_fn_33936T7�4
-�*
(�%
inputs��������� 
� "�����������:
__inference_loss_fn_0_34085'�

� 
� "� :
__inference_loss_fn_1_34096)�

� 
� "� �
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_31317�R�O
H�E
C�@
inputs4������������������������������������
� "H�E
>�;
04������������������������������������
� �
/__inference_max_pooling2d_3_layer_call_fn_31323�R�O
H�E
C�@
inputs4������������������������������������
� ";�84�������������������������������������
J__inference_max_pooling2d_4_layer_call_and_return_conditional_losses_31329�R�O
H�E
C�@
inputs4������������������������������������
� "H�E
>�;
04������������������������������������
� �
/__inference_max_pooling2d_4_layer_call_fn_31335�R�O
H�E
C�@
inputs4������������������������������������
� ";�84�������������������������������������
J__inference_max_pooling2d_5_layer_call_and_return_conditional_losses_31341�R�O
H�E
C�@
inputs4������������������������������������
� "H�E
>�;
04������������������������������������
� �
/__inference_max_pooling2d_5_layer_call_fn_31347�R�O
H�E
C�@
inputs4������������������������������������
� ";�84�������������������������������������
G__inference_sequential_3_layer_call_and_return_conditional_losses_30958����L�I
B�?
5�2
random_flip_1_input���������&
p

 
� "-�*
#� 
0���������&
� �
G__inference_sequential_3_layer_call_and_return_conditional_losses_30962}L�I
B�?
5�2
random_flip_1_input���������&
p 

 
� "-�*
#� 
0���������&
� �
G__inference_sequential_3_layer_call_and_return_conditional_losses_33791x���?�<
5�2
(�%
inputs���������&
p

 
� "-�*
#� 
0���������&
� �
G__inference_sequential_3_layer_call_and_return_conditional_losses_33795p?�<
5�2
(�%
inputs���������&
p 

 
� "-�*
#� 
0���������&
� �
,__inference_sequential_3_layer_call_fn_31302x���L�I
B�?
5�2
random_flip_1_input���������&
p

 
� " ����������&�
,__inference_sequential_3_layer_call_fn_31311pL�I
B�?
5�2
random_flip_1_input���������&
p 

 
� " ����������&�
,__inference_sequential_3_layer_call_fn_33806k���?�<
5�2
(�%
inputs���������&
p

 
� " ����������&�
,__inference_sequential_3_layer_call_fn_33811c?�<
5�2
(�%
inputs���������&
p 

 
� " ����������&�
G__inference_sequential_4_layer_call_and_return_conditional_losses_31702����!"#$%&'()*+,K�H
A�>
4�1
sequential_3_input���������&
p

 
� "%�"
�
0���������
� �
G__inference_sequential_4_layer_call_and_return_conditional_losses_31761�!"#$%&'()*+,K�H
A�>
4�1
sequential_3_input���������&
p 

 
� "%�"
�
0���������
� �
G__inference_sequential_4_layer_call_and_return_conditional_losses_33318|���!"#$%&'()*+,?�<
5�2
(�%
inputs���������&
p

 
� "%�"
�
0���������
� �
G__inference_sequential_4_layer_call_and_return_conditional_losses_33389v!"#$%&'()*+,?�<
5�2
(�%
inputs���������&
p 

 
� "%�"
�
0���������
� �
,__inference_sequential_4_layer_call_fn_31862{���!"#$%&'()*+,K�H
A�>
4�1
sequential_3_input���������&
p

 
� "�����������
,__inference_sequential_4_layer_call_fn_31950u!"#$%&'()*+,K�H
A�>
4�1
sequential_3_input���������&
p 

 
� "�����������
,__inference_sequential_4_layer_call_fn_33424o���!"#$%&'()*+,?�<
5�2
(�%
inputs���������&
p

 
� "�����������
,__inference_sequential_4_layer_call_fn_33453i!"#$%&'()*+,?�<
5�2
(�%
inputs���������&
p 

 
� "�����������
G__inference_sequential_5_layer_call_and_return_conditional_losses_32075����!"#$%&'()*+,K�H
A�>
4�1
sequential_4_input���������&
p

 
� "%�"
�
0���������
� �
G__inference_sequential_5_layer_call_and_return_conditional_losses_32117�!"#$%&'()*+,K�H
A�>
4�1
sequential_4_input���������&
p 

 
� "%�"
�
0���������
� �
G__inference_sequential_5_layer_call_and_return_conditional_losses_32739|���!"#$%&'()*+,?�<
5�2
(�%
inputs���������&
p

 
� "%�"
�
0���������
� �
G__inference_sequential_5_layer_call_and_return_conditional_losses_32811v!"#$%&'()*+,?�<
5�2
(�%
inputs���������&
p 

 
� "%�"
�
0���������
� �
,__inference_sequential_5_layer_call_fn_32201{���!"#$%&'()*+,K�H
A�>
4�1
sequential_4_input���������&
p

 
� "�����������
,__inference_sequential_5_layer_call_fn_32272u!"#$%&'()*+,K�H
A�>
4�1
sequential_4_input���������&
p 

 
� "�����������
,__inference_sequential_5_layer_call_fn_32846o���!"#$%&'()*+,?�<
5�2
(�%
inputs���������&
p

 
� "�����������
,__inference_sequential_5_layer_call_fn_32875i!"#$%&'()*+,?�<
5�2
(�%
inputs���������&
p 

 
� "�����������
#__inference_signature_wrapper_32315�!"#$%&'()*+,Y�V
� 
O�L
J
sequential_4_input4�1
sequential_4_input���������&"5�2
0
	softmax_1#� 
	softmax_1����������
D__inference_softmax_1_layer_call_and_return_conditional_losses_33458\3�0
)�&
 �
inputs���������

 
� "%�"
�
0���������
� |
)__inference_softmax_1_layer_call_fn_33463O3�0
)�&
 �
inputs���������

 
� "����������