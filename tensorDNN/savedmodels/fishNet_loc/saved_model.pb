ё╞&
Ї┼
D
AddV2
x"T
y"T
z"T"
Ttype:
2	АР
B
AssignVariableOp
resource
value"dtype"
dtypetypeИ
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
Ы
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
В
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
delete_old_dirsbool(И
?
Mul
x"T
y"T
z"T"
Ttype:
2	Р
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
dtypetypeИ
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
list(type)(0И
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0И
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
9
Softmax
logits"T
softmax"T"
Ttype:
2
╛
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
executor_typestring И
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
Ц
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 И"serve*2.6.02v2.6.0-rc2-32-g919f693420e8хй#
~
conv2d/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv2d/kernel
w
!conv2d/kernel/Read/ReadVariableOpReadVariableOpconv2d/kernel*&
_output_shapes
:*
dtype0
n
conv2d/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv2d/bias
g
conv2d/bias/Read/ReadVariableOpReadVariableOpconv2d/bias*
_output_shapes
:*
dtype0
В
conv2d_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameconv2d_1/kernel
{
#conv2d_1/kernel/Read/ReadVariableOpReadVariableOpconv2d_1/kernel*&
_output_shapes
:*
dtype0
r
conv2d_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv2d_1/bias
k
!conv2d_1/bias/Read/ReadVariableOpReadVariableOpconv2d_1/bias*
_output_shapes
:*
dtype0
В
conv2d_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape: * 
shared_nameconv2d_2/kernel
{
#conv2d_2/kernel/Read/ReadVariableOpReadVariableOpconv2d_2/kernel*&
_output_shapes
: *
dtype0
r
conv2d_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameconv2d_2/bias
k
!conv2d_2/bias/Read/ReadVariableOpReadVariableOpconv2d_2/bias*
_output_shapes
: *
dtype0
В
conv2d_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape: @* 
shared_nameconv2d_3/kernel
{
#conv2d_3/kernel/Read/ReadVariableOpReadVariableOpconv2d_3/kernel*&
_output_shapes
: @*
dtype0
r
conv2d_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameconv2d_3/bias
k
!conv2d_3/bias/Read/ReadVariableOpReadVariableOpconv2d_3/bias*
_output_shapes
:@*
dtype0
В
conv2d_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:@P* 
shared_nameconv2d_4/kernel
{
#conv2d_4/kernel/Read/ReadVariableOpReadVariableOpconv2d_4/kernel*&
_output_shapes
:@P*
dtype0
r
conv2d_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*
shared_nameconv2d_4/bias
k
!conv2d_4/bias/Read/ReadVariableOpReadVariableOpconv2d_4/bias*
_output_shapes
:P*
dtype0
v
dense/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
аю*
shared_namedense/kernel
o
 dense/kernel/Read/ReadVariableOpReadVariableOpdense/kernel* 
_output_shapes
:
аю*
dtype0
m

dense/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:ю*
shared_name
dense/bias
f
dense/bias/Read/ReadVariableOpReadVariableOp
dense/bias*
_output_shapes	
:ю*
dtype0
z
dense_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
юм*
shared_namedense_1/kernel
s
"dense_1/kernel/Read/ReadVariableOpReadVariableOpdense_1/kernel* 
_output_shapes
:
юм*
dtype0
q
dense_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:м*
shared_namedense_1/bias
j
 dense_1/bias/Read/ReadVariableOpReadVariableOpdense_1/bias*
_output_shapes	
:м*
dtype0
y
dense_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	мd*
shared_namedense_2/kernel
r
"dense_2/kernel/Read/ReadVariableOpReadVariableOpdense_2/kernel*
_output_shapes
:	мd*
dtype0
p
dense_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*
shared_namedense_2/bias
i
 dense_2/bias/Read/ReadVariableOpReadVariableOpdense_2/bias*
_output_shapes
:d*
dtype0
x
dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:d2*
shared_namedense_3/kernel
q
"dense_3/kernel/Read/ReadVariableOpReadVariableOpdense_3/kernel*
_output_shapes

:d2*
dtype0
p
dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_3/bias
i
 dense_3/bias/Read/ReadVariableOpReadVariableOpdense_3/bias*
_output_shapes
:2*
dtype0
v
output/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*
shared_nameoutput/kernel
o
!output/kernel/Read/ReadVariableOpReadVariableOpoutput/kernel*
_output_shapes

:2*
dtype0
n
output/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameoutput/bias
g
output/bias/Read/ReadVariableOpReadVariableOpoutput/bias*
_output_shapes
:*
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
М
Adam/conv2d/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv2d/kernel/m
Е
(Adam/conv2d/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv2d/kernel/m*&
_output_shapes
:*
dtype0
|
Adam/conv2d/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/conv2d/bias/m
u
&Adam/conv2d/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv2d/bias/m*
_output_shapes
:*
dtype0
Р
Adam/conv2d_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv2d_1/kernel/m
Й
*Adam/conv2d_1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_1/kernel/m*&
_output_shapes
:*
dtype0
А
Adam/conv2d_1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv2d_1/bias/m
y
(Adam/conv2d_1/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_1/bias/m*
_output_shapes
:*
dtype0
Р
Adam/conv2d_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: *'
shared_nameAdam/conv2d_2/kernel/m
Й
*Adam/conv2d_2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_2/kernel/m*&
_output_shapes
: *
dtype0
А
Adam/conv2d_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/conv2d_2/bias/m
y
(Adam/conv2d_2/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_2/bias/m*
_output_shapes
: *
dtype0
Р
Adam/conv2d_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape: @*'
shared_nameAdam/conv2d_3/kernel/m
Й
*Adam/conv2d_3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_3/kernel/m*&
_output_shapes
: @*
dtype0
А
Adam/conv2d_3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*%
shared_nameAdam/conv2d_3/bias/m
y
(Adam/conv2d_3/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_3/bias/m*
_output_shapes
:@*
dtype0
Р
Adam/conv2d_4/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:@P*'
shared_nameAdam/conv2d_4/kernel/m
Й
*Adam/conv2d_4/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_4/kernel/m*&
_output_shapes
:@P*
dtype0
А
Adam/conv2d_4/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*%
shared_nameAdam/conv2d_4/bias/m
y
(Adam/conv2d_4/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv2d_4/bias/m*
_output_shapes
:P*
dtype0
Д
Adam/dense/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
аю*$
shared_nameAdam/dense/kernel/m
}
'Adam/dense/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense/kernel/m* 
_output_shapes
:
аю*
dtype0
{
Adam/dense/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:ю*"
shared_nameAdam/dense/bias/m
t
%Adam/dense/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense/bias/m*
_output_shapes	
:ю*
dtype0
И
Adam/dense_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
юм*&
shared_nameAdam/dense_1/kernel/m
Б
)Adam/dense_1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_1/kernel/m* 
_output_shapes
:
юм*
dtype0

Adam/dense_1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:м*$
shared_nameAdam/dense_1/bias/m
x
'Adam/dense_1/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_1/bias/m*
_output_shapes	
:м*
dtype0
З
Adam/dense_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	мd*&
shared_nameAdam/dense_2/kernel/m
А
)Adam/dense_2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_2/kernel/m*
_output_shapes
:	мd*
dtype0
~
Adam/dense_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*$
shared_nameAdam/dense_2/bias/m
w
'Adam/dense_2/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_2/bias/m*
_output_shapes
:d*
dtype0
Ж
Adam/dense_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:d2*&
shared_nameAdam/dense_3/kernel/m

)Adam/dense_3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_3/kernel/m*
_output_shapes

:d2*
dtype0
~
Adam/dense_3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/dense_3/bias/m
w
'Adam/dense_3/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_3/bias/m*
_output_shapes
:2*
dtype0
Д
Adam/output/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*%
shared_nameAdam/output/kernel/m
}
(Adam/output/kernel/m/Read/ReadVariableOpReadVariableOpAdam/output/kernel/m*
_output_shapes

:2*
dtype0
|
Adam/output/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/output/bias/m
u
&Adam/output/bias/m/Read/ReadVariableOpReadVariableOpAdam/output/bias/m*
_output_shapes
:*
dtype0
М
Adam/conv2d/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv2d/kernel/v
Е
(Adam/conv2d/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv2d/kernel/v*&
_output_shapes
:*
dtype0
|
Adam/conv2d/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/conv2d/bias/v
u
&Adam/conv2d/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv2d/bias/v*
_output_shapes
:*
dtype0
Р
Adam/conv2d_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv2d_1/kernel/v
Й
*Adam/conv2d_1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_1/kernel/v*&
_output_shapes
:*
dtype0
А
Adam/conv2d_1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv2d_1/bias/v
y
(Adam/conv2d_1/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_1/bias/v*
_output_shapes
:*
dtype0
Р
Adam/conv2d_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: *'
shared_nameAdam/conv2d_2/kernel/v
Й
*Adam/conv2d_2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_2/kernel/v*&
_output_shapes
: *
dtype0
А
Adam/conv2d_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/conv2d_2/bias/v
y
(Adam/conv2d_2/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_2/bias/v*
_output_shapes
: *
dtype0
Р
Adam/conv2d_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape: @*'
shared_nameAdam/conv2d_3/kernel/v
Й
*Adam/conv2d_3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_3/kernel/v*&
_output_shapes
: @*
dtype0
А
Adam/conv2d_3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*%
shared_nameAdam/conv2d_3/bias/v
y
(Adam/conv2d_3/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_3/bias/v*
_output_shapes
:@*
dtype0
Р
Adam/conv2d_4/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:@P*'
shared_nameAdam/conv2d_4/kernel/v
Й
*Adam/conv2d_4/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_4/kernel/v*&
_output_shapes
:@P*
dtype0
А
Adam/conv2d_4/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:P*%
shared_nameAdam/conv2d_4/bias/v
y
(Adam/conv2d_4/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv2d_4/bias/v*
_output_shapes
:P*
dtype0
Д
Adam/dense/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
аю*$
shared_nameAdam/dense/kernel/v
}
'Adam/dense/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense/kernel/v* 
_output_shapes
:
аю*
dtype0
{
Adam/dense/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:ю*"
shared_nameAdam/dense/bias/v
t
%Adam/dense/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense/bias/v*
_output_shapes	
:ю*
dtype0
И
Adam/dense_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
юм*&
shared_nameAdam/dense_1/kernel/v
Б
)Adam/dense_1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_1/kernel/v* 
_output_shapes
:
юм*
dtype0

Adam/dense_1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:м*$
shared_nameAdam/dense_1/bias/v
x
'Adam/dense_1/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_1/bias/v*
_output_shapes	
:м*
dtype0
З
Adam/dense_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	мd*&
shared_nameAdam/dense_2/kernel/v
А
)Adam/dense_2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_2/kernel/v*
_output_shapes
:	мd*
dtype0
~
Adam/dense_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*$
shared_nameAdam/dense_2/bias/v
w
'Adam/dense_2/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_2/bias/v*
_output_shapes
:d*
dtype0
Ж
Adam/dense_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:d2*&
shared_nameAdam/dense_3/kernel/v

)Adam/dense_3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_3/kernel/v*
_output_shapes

:d2*
dtype0
~
Adam/dense_3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*$
shared_nameAdam/dense_3/bias/v
w
'Adam/dense_3/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_3/bias/v*
_output_shapes
:2*
dtype0
Д
Adam/output/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:2*%
shared_nameAdam/output/kernel/v
}
(Adam/output/kernel/v/Read/ReadVariableOpReadVariableOpAdam/output/kernel/v*
_output_shapes

:2*
dtype0
|
Adam/output/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/output/bias/v
u
&Adam/output/bias/v/Read/ReadVariableOpReadVariableOpAdam/output/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
╩В
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*ДВ
value∙БBїБ BэБ
ї
layer-0
layer-1
layer_with_weights-0
layer-2
layer-3
layer_with_weights-1
layer-4
layer-5
layer_with_weights-2
layer-6
layer-7
	layer_with_weights-3
	layer-8

layer_with_weights-4

layer-9
layer-10
layer-11
layer_with_weights-5
layer-12
layer_with_weights-6
layer-13
layer-14
layer_with_weights-7
layer-15
layer_with_weights-8
layer-16
layer-17
layer_with_weights-9
layer-18
	optimizer
	variables
regularization_losses
trainable_variables
	keras_api

signatures
y
layer-0
layer-1
layer-2
	variables
regularization_losses
trainable_variables
 	keras_api
R
!	variables
"regularization_losses
#trainable_variables
$	keras_api
h

%kernel
&bias
'	variables
(regularization_losses
)trainable_variables
*	keras_api
R
+	variables
,regularization_losses
-trainable_variables
.	keras_api
h

/kernel
0bias
1	variables
2regularization_losses
3trainable_variables
4	keras_api
R
5	variables
6regularization_losses
7trainable_variables
8	keras_api
h

9kernel
:bias
;	variables
<regularization_losses
=trainable_variables
>	keras_api
R
?	variables
@regularization_losses
Atrainable_variables
B	keras_api
h

Ckernel
Dbias
E	variables
Fregularization_losses
Gtrainable_variables
H	keras_api
h

Ikernel
Jbias
K	variables
Lregularization_losses
Mtrainable_variables
N	keras_api
R
O	variables
Pregularization_losses
Qtrainable_variables
R	keras_api
R
S	variables
Tregularization_losses
Utrainable_variables
V	keras_api
h

Wkernel
Xbias
Y	variables
Zregularization_losses
[trainable_variables
\	keras_api
h

]kernel
^bias
_	variables
`regularization_losses
atrainable_variables
b	keras_api
R
c	variables
dregularization_losses
etrainable_variables
f	keras_api
h

gkernel
hbias
i	variables
jregularization_losses
ktrainable_variables
l	keras_api
h

mkernel
nbias
o	variables
pregularization_losses
qtrainable_variables
r	keras_api
R
s	variables
tregularization_losses
utrainable_variables
v	keras_api
h

wkernel
xbias
y	variables
zregularization_losses
{trainable_variables
|	keras_api
╥
}iter

~beta_1

beta_2

Аdecay
Бlearning_rate%mТ&mУ/mФ0mХ9mЦ:mЧCmШDmЩImЪJmЫWmЬXmЭ]mЮ^mЯgmаhmбmmвnmгwmдxmе%vж&vз/vи0vй9vк:vлCvмDvнIvоJvпWv░Xv▒]v▓^v│gv┤hv╡mv╢nv╖wv╕xv╣
Ц
%0
&1
/2
03
94
:5
C6
D7
I8
J9
W10
X11
]12
^13
g14
h15
m16
n17
w18
x19
 
Ц
%0
&1
/2
03
94
:5
C6
D7
I8
J9
W10
X11
]12
^13
g14
h15
m16
n17
w18
x19
▓
Вmetrics
Гlayers
	variables
Дlayer_metrics
Еnon_trainable_variables
regularization_losses
 Жlayer_regularization_losses
trainable_variables
 
a
	З_rng
И	variables
Йregularization_losses
Кtrainable_variables
Л	keras_api
a
	М_rng
Н	variables
Оregularization_losses
Пtrainable_variables
Р	keras_api
a
	С_rng
Т	variables
Уregularization_losses
Фtrainable_variables
Х	keras_api
 
 
 
▓
Цmetrics
Чlayers
	variables
Шlayer_metrics
Щnon_trainable_variables
regularization_losses
 Ъlayer_regularization_losses
trainable_variables
 
 
 
▓
Ыmetrics
Ьlayers
Эlayer_metrics
!	variables
Юnon_trainable_variables
"regularization_losses
 Яlayer_regularization_losses
#trainable_variables
YW
VARIABLE_VALUEconv2d/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv2d/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE

%0
&1
 

%0
&1
▓
аmetrics
бlayers
вlayer_metrics
'	variables
гnon_trainable_variables
(regularization_losses
 дlayer_regularization_losses
)trainable_variables
 
 
 
▓
еmetrics
жlayers
зlayer_metrics
+	variables
иnon_trainable_variables
,regularization_losses
 йlayer_regularization_losses
-trainable_variables
[Y
VARIABLE_VALUEconv2d_1/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEconv2d_1/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE

/0
01
 

/0
01
▓
кmetrics
лlayers
мlayer_metrics
1	variables
нnon_trainable_variables
2regularization_losses
 оlayer_regularization_losses
3trainable_variables
 
 
 
▓
пmetrics
░layers
▒layer_metrics
5	variables
▓non_trainable_variables
6regularization_losses
 │layer_regularization_losses
7trainable_variables
[Y
VARIABLE_VALUEconv2d_2/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEconv2d_2/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE

90
:1
 

90
:1
▓
┤metrics
╡layers
╢layer_metrics
;	variables
╖non_trainable_variables
<regularization_losses
 ╕layer_regularization_losses
=trainable_variables
 
 
 
▓
╣metrics
║layers
╗layer_metrics
?	variables
╝non_trainable_variables
@regularization_losses
 ╜layer_regularization_losses
Atrainable_variables
[Y
VARIABLE_VALUEconv2d_3/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEconv2d_3/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE

C0
D1
 

C0
D1
▓
╛metrics
┐layers
└layer_metrics
E	variables
┴non_trainable_variables
Fregularization_losses
 ┬layer_regularization_losses
Gtrainable_variables
[Y
VARIABLE_VALUEconv2d_4/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEconv2d_4/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE

I0
J1
 

I0
J1
▓
├metrics
─layers
┼layer_metrics
K	variables
╞non_trainable_variables
Lregularization_losses
 ╟layer_regularization_losses
Mtrainable_variables
 
 
 
▓
╚metrics
╔layers
╩layer_metrics
O	variables
╦non_trainable_variables
Pregularization_losses
 ╠layer_regularization_losses
Qtrainable_variables
 
 
 
▓
═metrics
╬layers
╧layer_metrics
S	variables
╨non_trainable_variables
Tregularization_losses
 ╤layer_regularization_losses
Utrainable_variables
XV
VARIABLE_VALUEdense/kernel6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUE
TR
VARIABLE_VALUE
dense/bias4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUE

W0
X1
 

W0
X1
▓
╥metrics
╙layers
╘layer_metrics
Y	variables
╒non_trainable_variables
Zregularization_losses
 ╓layer_regularization_losses
[trainable_variables
ZX
VARIABLE_VALUEdense_1/kernel6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_1/bias4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUE

]0
^1
 

]0
^1
▓
╫metrics
╪layers
┘layer_metrics
_	variables
┌non_trainable_variables
`regularization_losses
 █layer_regularization_losses
atrainable_variables
 
 
 
▓
▄metrics
▌layers
▐layer_metrics
c	variables
▀non_trainable_variables
dregularization_losses
 рlayer_regularization_losses
etrainable_variables
ZX
VARIABLE_VALUEdense_2/kernel6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_2/bias4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUE

g0
h1
 

g0
h1
▓
сmetrics
тlayers
уlayer_metrics
i	variables
фnon_trainable_variables
jregularization_losses
 хlayer_regularization_losses
ktrainable_variables
ZX
VARIABLE_VALUEdense_3/kernel6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_3/bias4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUE

m0
n1
 

m0
n1
▓
цmetrics
чlayers
шlayer_metrics
o	variables
щnon_trainable_variables
pregularization_losses
 ъlayer_regularization_losses
qtrainable_variables
 
 
 
▓
ыmetrics
ьlayers
эlayer_metrics
s	variables
юnon_trainable_variables
tregularization_losses
 яlayer_regularization_losses
utrainable_variables
YW
VARIABLE_VALUEoutput/kernel6layer_with_weights-9/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEoutput/bias4layer_with_weights-9/bias/.ATTRIBUTES/VARIABLE_VALUE

w0
x1
 

w0
x1
▓
Ёmetrics
ёlayers
Єlayer_metrics
y	variables
єnon_trainable_variables
zregularization_losses
 Їlayer_regularization_losses
{trainable_variables
HF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
ZX
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE

ї0
Ў1
О
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
14
15
16
17
18
 
 
 

ў
_state_var
 
 
 
╡
°metrics
∙layers
·layer_metrics
И	variables
√non_trainable_variables
Йregularization_losses
 №layer_regularization_losses
Кtrainable_variables

¤
_state_var
 
 
 
╡
■metrics
 layers
Аlayer_metrics
Н	variables
Бnon_trainable_variables
Оregularization_losses
 Вlayer_regularization_losses
Пtrainable_variables

Г
_state_var
 
 
 
╡
Дmetrics
Еlayers
Жlayer_metrics
Т	variables
Зnon_trainable_variables
Уregularization_losses
 Иlayer_regularization_losses
Фtrainable_variables
 

0
1
2
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

Йtotal

Кcount
Л	variables
М	keras_api
I

Нtotal

Оcount
П
_fn_kwargs
Р	variables
С	keras_api
XV
VARIABLE_VALUEVariable:layer-0/layer-0/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUE
 
 
 
 
 
ZX
VARIABLE_VALUE
Variable_1:layer-0/layer-1/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUE
 
 
 
 
 
ZX
VARIABLE_VALUE
Variable_2:layer-0/layer-2/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUE
 
 
 
 
 
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

Й0
К1

Л	variables
QO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE
 

Н0
О1

Р	variables
|z
VARIABLE_VALUEAdam/conv2d/kernel/mRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv2d/bias/mPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/conv2d_1/kernel/mRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/conv2d_1/bias/mPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/conv2d_2/kernel/mRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/conv2d_2/bias/mPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/conv2d_3/kernel/mRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/conv2d_3/bias/mPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/conv2d_4/kernel/mRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/conv2d_4/bias/mPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
{y
VARIABLE_VALUEAdam/dense/kernel/mRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
wu
VARIABLE_VALUEAdam/dense/bias/mPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_1/kernel/mRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_1/bias/mPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_2/kernel/mRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_2/bias/mPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_3/kernel/mRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_3/bias/mPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/output/kernel/mRlayer_with_weights-9/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/output/bias/mPlayer_with_weights-9/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/conv2d/kernel/vRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv2d/bias/vPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/conv2d_1/kernel/vRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/conv2d_1/bias/vPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/conv2d_2/kernel/vRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/conv2d_2/bias/vPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/conv2d_3/kernel/vRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/conv2d_3/bias/vPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/conv2d_4/kernel/vRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUEAdam/conv2d_4/bias/vPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
{y
VARIABLE_VALUEAdam/dense/kernel/vRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
wu
VARIABLE_VALUEAdam/dense/bias/vPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_1/kernel/vRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_1/bias/vPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_2/kernel/vRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_2/bias/vPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_3/kernel/vRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_3/bias/vPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
|z
VARIABLE_VALUEAdam/output/kernel/vRlayer_with_weights-9/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/output/bias/vPlayer_with_weights-9/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
У
 serving_default_sequential_inputPlaceholder*/
_output_shapes
:         &*
dtype0*$
shape:         &
О
StatefulPartitionedCallStatefulPartitionedCall serving_default_sequential_inputconv2d/kernelconv2d/biasconv2d_1/kernelconv2d_1/biasconv2d_2/kernelconv2d_2/biasconv2d_3/kernelconv2d_3/biasconv2d_4/kernelconv2d_4/biasdense/kernel
dense/biasdense_1/kerneldense_1/biasdense_2/kerneldense_2/biasdense_3/kerneldense_3/biasoutput/kerneloutput/bias* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8В *,
f'R%
#__inference_signature_wrapper_77313
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
е
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename!conv2d/kernel/Read/ReadVariableOpconv2d/bias/Read/ReadVariableOp#conv2d_1/kernel/Read/ReadVariableOp!conv2d_1/bias/Read/ReadVariableOp#conv2d_2/kernel/Read/ReadVariableOp!conv2d_2/bias/Read/ReadVariableOp#conv2d_3/kernel/Read/ReadVariableOp!conv2d_3/bias/Read/ReadVariableOp#conv2d_4/kernel/Read/ReadVariableOp!conv2d_4/bias/Read/ReadVariableOp dense/kernel/Read/ReadVariableOpdense/bias/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOp dense_1/bias/Read/ReadVariableOp"dense_2/kernel/Read/ReadVariableOp dense_2/bias/Read/ReadVariableOp"dense_3/kernel/Read/ReadVariableOp dense_3/bias/Read/ReadVariableOp!output/kernel/Read/ReadVariableOpoutput/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOpVariable/Read/ReadVariableOpVariable_1/Read/ReadVariableOpVariable_2/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOp(Adam/conv2d/kernel/m/Read/ReadVariableOp&Adam/conv2d/bias/m/Read/ReadVariableOp*Adam/conv2d_1/kernel/m/Read/ReadVariableOp(Adam/conv2d_1/bias/m/Read/ReadVariableOp*Adam/conv2d_2/kernel/m/Read/ReadVariableOp(Adam/conv2d_2/bias/m/Read/ReadVariableOp*Adam/conv2d_3/kernel/m/Read/ReadVariableOp(Adam/conv2d_3/bias/m/Read/ReadVariableOp*Adam/conv2d_4/kernel/m/Read/ReadVariableOp(Adam/conv2d_4/bias/m/Read/ReadVariableOp'Adam/dense/kernel/m/Read/ReadVariableOp%Adam/dense/bias/m/Read/ReadVariableOp)Adam/dense_1/kernel/m/Read/ReadVariableOp'Adam/dense_1/bias/m/Read/ReadVariableOp)Adam/dense_2/kernel/m/Read/ReadVariableOp'Adam/dense_2/bias/m/Read/ReadVariableOp)Adam/dense_3/kernel/m/Read/ReadVariableOp'Adam/dense_3/bias/m/Read/ReadVariableOp(Adam/output/kernel/m/Read/ReadVariableOp&Adam/output/bias/m/Read/ReadVariableOp(Adam/conv2d/kernel/v/Read/ReadVariableOp&Adam/conv2d/bias/v/Read/ReadVariableOp*Adam/conv2d_1/kernel/v/Read/ReadVariableOp(Adam/conv2d_1/bias/v/Read/ReadVariableOp*Adam/conv2d_2/kernel/v/Read/ReadVariableOp(Adam/conv2d_2/bias/v/Read/ReadVariableOp*Adam/conv2d_3/kernel/v/Read/ReadVariableOp(Adam/conv2d_3/bias/v/Read/ReadVariableOp*Adam/conv2d_4/kernel/v/Read/ReadVariableOp(Adam/conv2d_4/bias/v/Read/ReadVariableOp'Adam/dense/kernel/v/Read/ReadVariableOp%Adam/dense/bias/v/Read/ReadVariableOp)Adam/dense_1/kernel/v/Read/ReadVariableOp'Adam/dense_1/bias/v/Read/ReadVariableOp)Adam/dense_2/kernel/v/Read/ReadVariableOp'Adam/dense_2/bias/v/Read/ReadVariableOp)Adam/dense_3/kernel/v/Read/ReadVariableOp'Adam/dense_3/bias/v/Read/ReadVariableOp(Adam/output/kernel/v/Read/ReadVariableOp&Adam/output/bias/v/Read/ReadVariableOpConst*U
TinN
L2J				*
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
GPU 2J 8В *'
f"R 
__inference__traced_save_79277
А
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameconv2d/kernelconv2d/biasconv2d_1/kernelconv2d_1/biasconv2d_2/kernelconv2d_2/biasconv2d_3/kernelconv2d_3/biasconv2d_4/kernelconv2d_4/biasdense/kernel
dense/biasdense_1/kerneldense_1/biasdense_2/kerneldense_2/biasdense_3/kerneldense_3/biasoutput/kerneloutput/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_rateVariable
Variable_1
Variable_2totalcounttotal_1count_1Adam/conv2d/kernel/mAdam/conv2d/bias/mAdam/conv2d_1/kernel/mAdam/conv2d_1/bias/mAdam/conv2d_2/kernel/mAdam/conv2d_2/bias/mAdam/conv2d_3/kernel/mAdam/conv2d_3/bias/mAdam/conv2d_4/kernel/mAdam/conv2d_4/bias/mAdam/dense/kernel/mAdam/dense/bias/mAdam/dense_1/kernel/mAdam/dense_1/bias/mAdam/dense_2/kernel/mAdam/dense_2/bias/mAdam/dense_3/kernel/mAdam/dense_3/bias/mAdam/output/kernel/mAdam/output/bias/mAdam/conv2d/kernel/vAdam/conv2d/bias/vAdam/conv2d_1/kernel/vAdam/conv2d_1/bias/vAdam/conv2d_2/kernel/vAdam/conv2d_2/bias/vAdam/conv2d_3/kernel/vAdam/conv2d_3/bias/vAdam/conv2d_4/kernel/vAdam/conv2d_4/bias/vAdam/dense/kernel/vAdam/dense/bias/vAdam/dense_1/kernel/vAdam/dense_1/bias/vAdam/dense_2/kernel/vAdam/dense_2/bias/vAdam/dense_3/kernel/vAdam/dense_3/bias/vAdam/output/kernel/vAdam/output/bias/v*T
TinM
K2I*
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
GPU 2J 8В **
f%R#
!__inference__traced_restore_79503ўЁ 
▀}
Ї
G__inference_sequential_1_layer_call_and_return_conditional_losses_77236
sequential_input
sequential_77146:	
sequential_77148:	
sequential_77150:	&
conv2d_77154:
conv2d_77156:(
conv2d_1_77160:
conv2d_1_77162:(
conv2d_2_77166: 
conv2d_2_77168: (
conv2d_3_77172: @
conv2d_3_77174:@(
conv2d_4_77177:@P
conv2d_4_77179:P
dense_77184:
аю
dense_77186:	ю!
dense_1_77189:
юм
dense_1_77191:	м 
dense_2_77195:	мd
dense_2_77197:d
dense_3_77200:d2
dense_3_77202:2
output_77206:2
output_77208:
identityИвconv2d/StatefulPartitionedCallв conv2d_1/StatefulPartitionedCallв conv2d_2/StatefulPartitionedCallв conv2d_3/StatefulPartitionedCallв conv2d_4/StatefulPartitionedCallвdense/StatefulPartitionedCallв.dense/kernel/Regularizer/Square/ReadVariableOpвdense_1/StatefulPartitionedCallв0dense_1/kernel/Regularizer/Square/ReadVariableOpвdense_2/StatefulPartitionedCallв0dense_2/kernel/Regularizer/Square/ReadVariableOpвdense_3/StatefulPartitionedCallв0dense_3/kernel/Regularizer/Square/ReadVariableOpвdropout/StatefulPartitionedCallв!dropout_1/StatefulPartitionedCallвoutput/StatefulPartitionedCallв"sequential/StatefulPartitionedCall╝
"sequential/StatefulPartitionedCallStatefulPartitionedCallsequential_inputsequential_77146sequential_77148sequential_77150*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *N
fIRG
E__inference_sequential_layer_call_and_return_conditional_losses_761842$
"sequential/StatefulPartitionedCallЙ
rescaling_1/PartitionedCallPartitionedCall+sequential/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_rescaling_1_layer_call_and_return_conditional_losses_763282
rescaling_1/PartitionedCallн
conv2d/StatefulPartitionedCallStatefulPartitionedCall$rescaling_1/PartitionedCall:output:0conv2d_77154conv2d_77156*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_conv2d_layer_call_and_return_conditional_losses_763412 
conv2d/StatefulPartitionedCallЛ
max_pooling2d/PartitionedCallPartitionedCall'conv2d/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *Q
fLRJ
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_763512
max_pooling2d/PartitionedCall╣
 conv2d_1/StatefulPartitionedCallStatefulPartitionedCall&max_pooling2d/PartitionedCall:output:0conv2d_1_77160conv2d_1_77162*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_1_layer_call_and_return_conditional_losses_763642"
 conv2d_1/StatefulPartitionedCallУ
max_pooling2d_1/PartitionedCallPartitionedCall)conv2d_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_763742!
max_pooling2d_1/PartitionedCall╗
 conv2d_2/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_1/PartitionedCall:output:0conv2d_2_77166conv2d_2_77168*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_2_layer_call_and_return_conditional_losses_763872"
 conv2d_2/StatefulPartitionedCallУ
max_pooling2d_2/PartitionedCallPartitionedCall)conv2d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:          * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_763972!
max_pooling2d_2/PartitionedCall╗
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_2/PartitionedCall:output:0conv2d_3_77172conv2d_3_77174*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         @*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_764102"
 conv2d_3/StatefulPartitionedCall╝
 conv2d_4/StatefulPartitionedCallStatefulPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0conv2d_4_77177conv2d_4_77179*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_764272"
 conv2d_4/StatefulPartitionedCallУ
max_pooling2d_3/PartitionedCallPartitionedCall)conv2d_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_764372!
max_pooling2d_3/PartitionedCallє
flatten/PartitionedCallPartitionedCall(max_pooling2d_3/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         а* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_764452
flatten/PartitionedCallЭ
dense/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_77184dense_77186*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ю*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *I
fDRB
@__inference_dense_layer_call_and_return_conditional_losses_764642
dense/StatefulPartitionedCallн
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_77189dense_1_77191*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_1_layer_call_and_return_conditional_losses_764872!
dense_1/StatefulPartitionedCallЛ
dropout/StatefulPartitionedCallStatefulPartitionedCall(dense_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dropout_layer_call_and_return_conditional_losses_767112!
dropout/StatefulPartitionedCallо
dense_2/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0dense_2_77195dense_2_77197*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         d*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_765172!
dense_2/StatefulPartitionedCallо
dense_3/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0dense_3_77200dense_3_77202*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_765402!
dense_3/StatefulPartitionedCall▓
!dropout_1/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0 ^dropout/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *M
fHRF
D__inference_dropout_1_layer_call_and_return_conditional_losses_766682#
!dropout_1/StatefulPartitionedCallл
output/StatefulPartitionedCallStatefulPartitionedCall*dropout_1/StatefulPartitionedCall:output:0output_77206output_77208*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_output_layer_call_and_return_conditional_losses_765642 
output/StatefulPartitionedCallо
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_77184* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul┤
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_77189* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul│
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_77195*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul▓
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_77200*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mulВ
IdentityIdentity'output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identity╫
NoOpNoOp^conv2d/StatefulPartitionedCall!^conv2d_1/StatefulPartitionedCall!^conv2d_2/StatefulPartitionedCall!^conv2d_3/StatefulPartitionedCall!^conv2d_4/StatefulPartitionedCall^dense/StatefulPartitionedCall/^dense/kernel/Regularizer/Square/ReadVariableOp ^dense_1/StatefulPartitionedCall1^dense_1/kernel/Regularizer/Square/ReadVariableOp ^dense_2/StatefulPartitionedCall1^dense_2/kernel/Regularizer/Square/ReadVariableOp ^dense_3/StatefulPartitionedCall1^dense_3/kernel/Regularizer/Square/ReadVariableOp ^dropout/StatefulPartitionedCall"^dropout_1/StatefulPartitionedCall^output/StatefulPartitionedCall#^sequential/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*\
_input_shapesK
I:         &: : : : : : : : : : : : : : : : : : : : : : : 2@
conv2d/StatefulPartitionedCallconv2d/StatefulPartitionedCall2D
 conv2d_1/StatefulPartitionedCall conv2d_1/StatefulPartitionedCall2D
 conv2d_2/StatefulPartitionedCall conv2d_2/StatefulPartitionedCall2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2D
 conv2d_4/StatefulPartitionedCall conv2d_4/StatefulPartitionedCall2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2F
!dropout_1/StatefulPartitionedCall!dropout_1/StatefulPartitionedCall2@
output/StatefulPartitionedCalloutput/StatefulPartitionedCall2H
"sequential/StatefulPartitionedCall"sequential/StatefulPartitionedCall:a ]
/
_output_shapes
:         &
*
_user_specified_namesequential_input
З
Є
A__inference_output_layer_call_and_return_conditional_losses_78637

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpН
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2	
BiasAdda
SoftmaxSoftmaxBiasAdd:output:0*
T0*'
_output_shapes
:         2	
Softmaxl
IdentityIdentitySoftmax:softmax:0^NoOp*
T0*'
_output_shapes
:         2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
к
f
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_76277

inputs
identityн
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4                                    *
ksize
*
paddingVALID*
strides
2	
MaxPoolЗ
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
к
f
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_78319

inputs
identityн
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4                                    *
ksize
*
paddingVALID*
strides
2	
MaxPoolЗ
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
ю

/__inference_random_rotation_layer_call_fn_78767

inputs
unknown:	
identityИвStatefulPartitionedCallЄ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_random_rotation_layer_call_and_return_conditional_losses_760882
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         &2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ю
f
J__inference_random_rotation_layer_call_and_return_conditional_losses_78771

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ё
Х
%__inference_dense_layer_call_fn_78450

inputs
unknown:
аю
	unknown_0:	ю
identityИвStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ю*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *I
fDRB
@__inference_dense_layer_call_and_return_conditional_losses_764642
StatefulPartitionedCall|
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:         ю2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         а: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:         а
 
_user_specified_nameinputs
╝
C
'__inference_dropout_layer_call_fn_78504

inputs
identity┴
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dropout_layer_call_and_return_conditional_losses_764982
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:         м2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:         м:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
ш
K
/__inference_max_pooling2d_2_layer_call_fn_78354

inputs
identity╨
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:          * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_763972
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:          2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         	 :W S
/
_output_shapes
:         	 
 
_user_specified_nameinputs
щ
№
C__inference_conv2d_1_layer_call_and_return_conditional_losses_78304

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         *
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:         
 
_user_specified_nameinputs
╜
f
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_76437

inputs
identityТ
MaxPoolMaxPoolinputs*/
_output_shapes
:         P*
ksize
*
paddingVALID*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:         P2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         P:W S
/
_output_shapes
:         P
 
_user_specified_nameinputs
Ю
b
)__inference_dropout_1_layer_call_fn_78600

inputs
identityИвStatefulPartitionedCall┌
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *M
fHRF
D__inference_dropout_1_layer_call_and_return_conditional_losses_766682
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         22

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:         222
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
щ
ж
B__inference_dense_3_layer_call_and_return_conditional_losses_78590

inputs0
matmul_readvariableop_resource:d2-
biasadd_readvariableop_resource:2
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpв0dense_3/kernel/Regularizer/Square/ReadVariableOpН
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:d2*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:         22
Relu├
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mulm
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:         22

Identity▓
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         d: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:         d
 
_user_specified_nameinputs
▒
a
B__inference_dropout_layer_call_and_return_conditional_losses_76711

inputs
identityИc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n█╢?2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:         м2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape╡
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:         м*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЪЩЩ>2
dropout/GreaterEqual/y┐
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:         м2
dropout/GreaterEqualА
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:         м2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:         м2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:         м2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:         м:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
╜
f
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_78424

inputs
identityТ
MaxPoolMaxPoolinputs*/
_output_shapes
:         P*
ksize
*
paddingVALID*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:         P2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         P:W S
/
_output_shapes
:         P
 
_user_specified_nameinputs
╗
d
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_76351

inputs
identityТ
MaxPoolMaxPoolinputs*/
_output_shapes
:         *
ksize
*
paddingVALID*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:         2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ю
b
F__inference_rescaling_1_layer_call_and_return_conditional_losses_78244

inputs
identityU
Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *БАА;2
Cast/xY
Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2

Cast_1/xd
mulMulinputsCast/x:output:0*
T0*/
_output_shapes
:         &2
muli
addAddV2mul:z:0Cast_1/x:output:0*
T0*/
_output_shapes
:         &2
addc
IdentityIdentityadd:z:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
З
Є
A__inference_output_layer_call_and_return_conditional_losses_76564

inputs0
matmul_readvariableop_resource:2-
biasadd_readvariableop_resource:
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpН
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:2*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2	
BiasAdda
SoftmaxSoftmaxBiasAdd:output:0*
T0*'
_output_shapes
:         2	
Softmaxl
IdentityIdentitySoftmax:softmax:0^NoOp*
T0*'
_output_shapes
:         2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
Ы
╠
E__inference_sequential_layer_call_and_return_conditional_losses_76224
random_flip_input
random_flip_76214:	#
random_rotation_76217:	
random_zoom_76220:	
identityИв#random_flip/StatefulPartitionedCallв'random_rotation/StatefulPartitionedCallв#random_zoom/StatefulPartitionedCallЫ
#random_flip/StatefulPartitionedCallStatefulPartitionedCallrandom_flip_inputrandom_flip_76214*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_flip_layer_call_and_return_conditional_losses_761592%
#random_flip/StatefulPartitionedCall╞
'random_rotation/StatefulPartitionedCallStatefulPartitionedCall,random_flip/StatefulPartitionedCall:output:0random_rotation_76217*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_random_rotation_layer_call_and_return_conditional_losses_760882)
'random_rotation/StatefulPartitionedCall║
#random_zoom/StatefulPartitionedCallStatefulPartitionedCall0random_rotation/StatefulPartitionedCall:output:0random_zoom_76220*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_zoom_layer_call_and_return_conditional_losses_759532%
#random_zoom/StatefulPartitionedCallП
IdentityIdentity,random_zoom/StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         &2

Identity─
NoOpNoOp$^random_flip/StatefulPartitionedCall(^random_rotation/StatefulPartitionedCall$^random_zoom/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*4
_input_shapes#
!:         &: : : 2J
#random_flip/StatefulPartitionedCall#random_flip/StatefulPartitionedCall2R
'random_rotation/StatefulPartitionedCall'random_rotation/StatefulPartitionedCall2J
#random_zoom/StatefulPartitionedCall#random_zoom/StatefulPartitionedCall:b ^
/
_output_shapes
:         &
+
_user_specified_namerandom_flip_input
▒
a
B__inference_dropout_layer_call_and_return_conditional_losses_78526

inputs
identityИc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n█╢?2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:         м2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape╡
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:         м*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЪЩЩ>2
dropout/GreaterEqual/y┐
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:         м2
dropout/GreaterEqualА
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:         м2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:         м2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:         м2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:         м:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
╒
K
/__inference_max_pooling2d_1_layer_call_fn_78309

inputs
identityы
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4                                    * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_762552
PartitionedCallП
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
∙
й
B__inference_dense_1_layer_call_and_return_conditional_losses_76487

inputs2
matmul_readvariableop_resource:
юм.
biasadd_readvariableop_resource:	м
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpв0dense_1/kernel/Regularizer/Square/ReadVariableOpП
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:м*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         м2
Relu┼
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/muln
IdentityIdentityRelu:activations:0^NoOp*
T0*(
_output_shapes
:         м2

Identity▓
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_1/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         ю: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:         ю
 
_user_specified_nameinputs
╖
Е
,__inference_sequential_1_layer_call_fn_77409

inputs
unknown:	
	unknown_0:	
	unknown_1:	#
	unknown_2:
	unknown_3:#
	unknown_4:
	unknown_5:#
	unknown_6: 
	unknown_7: #
	unknown_8: @
	unknown_9:@$

unknown_10:@P

unknown_11:P

unknown_12:
аю

unknown_13:	ю

unknown_14:
юм

unknown_15:	м

unknown_16:	мd

unknown_17:d

unknown_18:d2

unknown_19:2

unknown_20:2

unknown_21:
identityИвStatefulPartitionedCallС
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21*#
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8В *P
fKRI
G__inference_sequential_1_layer_call_and_return_conditional_losses_769562
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*\
_input_shapesK
I:         &: : : : : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Тx
й

G__inference_sequential_1_layer_call_and_return_conditional_losses_77143
sequential_input&
conv2d_77061:
conv2d_77063:(
conv2d_1_77067:
conv2d_1_77069:(
conv2d_2_77073: 
conv2d_2_77075: (
conv2d_3_77079: @
conv2d_3_77081:@(
conv2d_4_77084:@P
conv2d_4_77086:P
dense_77091:
аю
dense_77093:	ю!
dense_1_77096:
юм
dense_1_77098:	м 
dense_2_77102:	мd
dense_2_77104:d
dense_3_77107:d2
dense_3_77109:2
output_77113:2
output_77115:
identityИвconv2d/StatefulPartitionedCallв conv2d_1/StatefulPartitionedCallв conv2d_2/StatefulPartitionedCallв conv2d_3/StatefulPartitionedCallв conv2d_4/StatefulPartitionedCallвdense/StatefulPartitionedCallв.dense/kernel/Regularizer/Square/ReadVariableOpвdense_1/StatefulPartitionedCallв0dense_1/kernel/Regularizer/Square/ReadVariableOpвdense_2/StatefulPartitionedCallв0dense_2/kernel/Regularizer/Square/ReadVariableOpвdense_3/StatefulPartitionedCallв0dense_3/kernel/Regularizer/Square/ReadVariableOpвoutput/StatefulPartitionedCallы
sequential/PartitionedCallPartitionedCallsequential_input*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *N
fIRG
E__inference_sequential_layer_call_and_return_conditional_losses_758282
sequential/PartitionedCallБ
rescaling_1/PartitionedCallPartitionedCall#sequential/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_rescaling_1_layer_call_and_return_conditional_losses_763282
rescaling_1/PartitionedCallн
conv2d/StatefulPartitionedCallStatefulPartitionedCall$rescaling_1/PartitionedCall:output:0conv2d_77061conv2d_77063*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_conv2d_layer_call_and_return_conditional_losses_763412 
conv2d/StatefulPartitionedCallЛ
max_pooling2d/PartitionedCallPartitionedCall'conv2d/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *Q
fLRJ
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_763512
max_pooling2d/PartitionedCall╣
 conv2d_1/StatefulPartitionedCallStatefulPartitionedCall&max_pooling2d/PartitionedCall:output:0conv2d_1_77067conv2d_1_77069*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_1_layer_call_and_return_conditional_losses_763642"
 conv2d_1/StatefulPartitionedCallУ
max_pooling2d_1/PartitionedCallPartitionedCall)conv2d_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_763742!
max_pooling2d_1/PartitionedCall╗
 conv2d_2/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_1/PartitionedCall:output:0conv2d_2_77073conv2d_2_77075*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_2_layer_call_and_return_conditional_losses_763872"
 conv2d_2/StatefulPartitionedCallУ
max_pooling2d_2/PartitionedCallPartitionedCall)conv2d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:          * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_763972!
max_pooling2d_2/PartitionedCall╗
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_2/PartitionedCall:output:0conv2d_3_77079conv2d_3_77081*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         @*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_764102"
 conv2d_3/StatefulPartitionedCall╝
 conv2d_4/StatefulPartitionedCallStatefulPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0conv2d_4_77084conv2d_4_77086*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_764272"
 conv2d_4/StatefulPartitionedCallУ
max_pooling2d_3/PartitionedCallPartitionedCall)conv2d_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_764372!
max_pooling2d_3/PartitionedCallє
flatten/PartitionedCallPartitionedCall(max_pooling2d_3/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         а* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_764452
flatten/PartitionedCallЭ
dense/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_77091dense_77093*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ю*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *I
fDRB
@__inference_dense_layer_call_and_return_conditional_losses_764642
dense/StatefulPartitionedCallн
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_77096dense_1_77098*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_1_layer_call_and_return_conditional_losses_764872!
dense_1/StatefulPartitionedCallє
dropout/PartitionedCallPartitionedCall(dense_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dropout_layer_call_and_return_conditional_losses_764982
dropout/PartitionedCallж
dense_2/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0dense_2_77102dense_2_77104*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         d*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_765172!
dense_2/StatefulPartitionedCallо
dense_3/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0dense_3_77107dense_3_77109*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_765402!
dense_3/StatefulPartitionedCall°
dropout_1/PartitionedCallPartitionedCall(dense_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *M
fHRF
D__inference_dropout_1_layer_call_and_return_conditional_losses_765512
dropout_1/PartitionedCallг
output/StatefulPartitionedCallStatefulPartitionedCall"dropout_1/PartitionedCall:output:0output_77113output_77115*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_output_layer_call_and_return_conditional_losses_765642 
output/StatefulPartitionedCallо
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_77091* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul┤
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_77096* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul│
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_77102*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul▓
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_77107*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mulВ
IdentityIdentity'output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identityь
NoOpNoOp^conv2d/StatefulPartitionedCall!^conv2d_1/StatefulPartitionedCall!^conv2d_2/StatefulPartitionedCall!^conv2d_3/StatefulPartitionedCall!^conv2d_4/StatefulPartitionedCall^dense/StatefulPartitionedCall/^dense/kernel/Regularizer/Square/ReadVariableOp ^dense_1/StatefulPartitionedCall1^dense_1/kernel/Regularizer/Square/ReadVariableOp ^dense_2/StatefulPartitionedCall1^dense_2/kernel/Regularizer/Square/ReadVariableOp ^dense_3/StatefulPartitionedCall1^dense_3/kernel/Regularizer/Square/ReadVariableOp^output/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*V
_input_shapesE
C:         &: : : : : : : : : : : : : : : : : : : : 2@
conv2d/StatefulPartitionedCallconv2d/StatefulPartitionedCall2D
 conv2d_1/StatefulPartitionedCall conv2d_1/StatefulPartitionedCall2D
 conv2d_2/StatefulPartitionedCall conv2d_2/StatefulPartitionedCall2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2D
 conv2d_4/StatefulPartitionedCall conv2d_4/StatefulPartitionedCall2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2@
output/StatefulPartitionedCalloutput/StatefulPartitionedCall:a ]
/
_output_shapes
:         &
*
_user_specified_namesequential_input
┴}
ъ
G__inference_sequential_1_layer_call_and_return_conditional_losses_76956

inputs
sequential_76866:	
sequential_76868:	
sequential_76870:	&
conv2d_76874:
conv2d_76876:(
conv2d_1_76880:
conv2d_1_76882:(
conv2d_2_76886: 
conv2d_2_76888: (
conv2d_3_76892: @
conv2d_3_76894:@(
conv2d_4_76897:@P
conv2d_4_76899:P
dense_76904:
аю
dense_76906:	ю!
dense_1_76909:
юм
dense_1_76911:	м 
dense_2_76915:	мd
dense_2_76917:d
dense_3_76920:d2
dense_3_76922:2
output_76926:2
output_76928:
identityИвconv2d/StatefulPartitionedCallв conv2d_1/StatefulPartitionedCallв conv2d_2/StatefulPartitionedCallв conv2d_3/StatefulPartitionedCallв conv2d_4/StatefulPartitionedCallвdense/StatefulPartitionedCallв.dense/kernel/Regularizer/Square/ReadVariableOpвdense_1/StatefulPartitionedCallв0dense_1/kernel/Regularizer/Square/ReadVariableOpвdense_2/StatefulPartitionedCallв0dense_2/kernel/Regularizer/Square/ReadVariableOpвdense_3/StatefulPartitionedCallв0dense_3/kernel/Regularizer/Square/ReadVariableOpвdropout/StatefulPartitionedCallв!dropout_1/StatefulPartitionedCallвoutput/StatefulPartitionedCallв"sequential/StatefulPartitionedCall▓
"sequential/StatefulPartitionedCallStatefulPartitionedCallinputssequential_76866sequential_76868sequential_76870*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *N
fIRG
E__inference_sequential_layer_call_and_return_conditional_losses_761842$
"sequential/StatefulPartitionedCallЙ
rescaling_1/PartitionedCallPartitionedCall+sequential/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_rescaling_1_layer_call_and_return_conditional_losses_763282
rescaling_1/PartitionedCallн
conv2d/StatefulPartitionedCallStatefulPartitionedCall$rescaling_1/PartitionedCall:output:0conv2d_76874conv2d_76876*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_conv2d_layer_call_and_return_conditional_losses_763412 
conv2d/StatefulPartitionedCallЛ
max_pooling2d/PartitionedCallPartitionedCall'conv2d/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *Q
fLRJ
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_763512
max_pooling2d/PartitionedCall╣
 conv2d_1/StatefulPartitionedCallStatefulPartitionedCall&max_pooling2d/PartitionedCall:output:0conv2d_1_76880conv2d_1_76882*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_1_layer_call_and_return_conditional_losses_763642"
 conv2d_1/StatefulPartitionedCallУ
max_pooling2d_1/PartitionedCallPartitionedCall)conv2d_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_763742!
max_pooling2d_1/PartitionedCall╗
 conv2d_2/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_1/PartitionedCall:output:0conv2d_2_76886conv2d_2_76888*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_2_layer_call_and_return_conditional_losses_763872"
 conv2d_2/StatefulPartitionedCallУ
max_pooling2d_2/PartitionedCallPartitionedCall)conv2d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:          * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_763972!
max_pooling2d_2/PartitionedCall╗
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_2/PartitionedCall:output:0conv2d_3_76892conv2d_3_76894*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         @*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_764102"
 conv2d_3/StatefulPartitionedCall╝
 conv2d_4/StatefulPartitionedCallStatefulPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0conv2d_4_76897conv2d_4_76899*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_764272"
 conv2d_4/StatefulPartitionedCallУ
max_pooling2d_3/PartitionedCallPartitionedCall)conv2d_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_764372!
max_pooling2d_3/PartitionedCallє
flatten/PartitionedCallPartitionedCall(max_pooling2d_3/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         а* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_764452
flatten/PartitionedCallЭ
dense/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_76904dense_76906*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ю*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *I
fDRB
@__inference_dense_layer_call_and_return_conditional_losses_764642
dense/StatefulPartitionedCallн
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_76909dense_1_76911*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_1_layer_call_and_return_conditional_losses_764872!
dense_1/StatefulPartitionedCallЛ
dropout/StatefulPartitionedCallStatefulPartitionedCall(dense_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dropout_layer_call_and_return_conditional_losses_767112!
dropout/StatefulPartitionedCallо
dense_2/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0dense_2_76915dense_2_76917*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         d*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_765172!
dense_2/StatefulPartitionedCallо
dense_3/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0dense_3_76920dense_3_76922*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_765402!
dense_3/StatefulPartitionedCall▓
!dropout_1/StatefulPartitionedCallStatefulPartitionedCall(dense_3/StatefulPartitionedCall:output:0 ^dropout/StatefulPartitionedCall*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *M
fHRF
D__inference_dropout_1_layer_call_and_return_conditional_losses_766682#
!dropout_1/StatefulPartitionedCallл
output/StatefulPartitionedCallStatefulPartitionedCall*dropout_1/StatefulPartitionedCall:output:0output_76926output_76928*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_output_layer_call_and_return_conditional_losses_765642 
output/StatefulPartitionedCallо
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_76904* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul┤
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_76909* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul│
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_76915*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul▓
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_76920*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mulВ
IdentityIdentity'output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identity╫
NoOpNoOp^conv2d/StatefulPartitionedCall!^conv2d_1/StatefulPartitionedCall!^conv2d_2/StatefulPartitionedCall!^conv2d_3/StatefulPartitionedCall!^conv2d_4/StatefulPartitionedCall^dense/StatefulPartitionedCall/^dense/kernel/Regularizer/Square/ReadVariableOp ^dense_1/StatefulPartitionedCall1^dense_1/kernel/Regularizer/Square/ReadVariableOp ^dense_2/StatefulPartitionedCall1^dense_2/kernel/Regularizer/Square/ReadVariableOp ^dense_3/StatefulPartitionedCall1^dense_3/kernel/Regularizer/Square/ReadVariableOp ^dropout/StatefulPartitionedCall"^dropout_1/StatefulPartitionedCall^output/StatefulPartitionedCall#^sequential/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*\
_input_shapesK
I:         &: : : : : : : : : : : : : : : : : : : : : : : 2@
conv2d/StatefulPartitionedCallconv2d/StatefulPartitionedCall2D
 conv2d_1/StatefulPartitionedCall conv2d_1/StatefulPartitionedCall2D
 conv2d_2/StatefulPartitionedCall conv2d_2/StatefulPartitionedCall2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2D
 conv2d_4/StatefulPartitionedCall conv2d_4/StatefulPartitionedCall2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2F
!dropout_1/StatefulPartitionedCall!dropout_1/StatefulPartitionedCall2@
output/StatefulPartitionedCalloutput/StatefulPartitionedCall2H
"sequential/StatefulPartitionedCall"sequential/StatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
■О
Ї
 __inference__wrapped_model_75802
sequential_inputL
2sequential_1_conv2d_conv2d_readvariableop_resource:A
3sequential_1_conv2d_biasadd_readvariableop_resource:N
4sequential_1_conv2d_1_conv2d_readvariableop_resource:C
5sequential_1_conv2d_1_biasadd_readvariableop_resource:N
4sequential_1_conv2d_2_conv2d_readvariableop_resource: C
5sequential_1_conv2d_2_biasadd_readvariableop_resource: N
4sequential_1_conv2d_3_conv2d_readvariableop_resource: @C
5sequential_1_conv2d_3_biasadd_readvariableop_resource:@N
4sequential_1_conv2d_4_conv2d_readvariableop_resource:@PC
5sequential_1_conv2d_4_biasadd_readvariableop_resource:PE
1sequential_1_dense_matmul_readvariableop_resource:
аюA
2sequential_1_dense_biasadd_readvariableop_resource:	юG
3sequential_1_dense_1_matmul_readvariableop_resource:
юмC
4sequential_1_dense_1_biasadd_readvariableop_resource:	мF
3sequential_1_dense_2_matmul_readvariableop_resource:	мdB
4sequential_1_dense_2_biasadd_readvariableop_resource:dE
3sequential_1_dense_3_matmul_readvariableop_resource:d2B
4sequential_1_dense_3_biasadd_readvariableop_resource:2D
2sequential_1_output_matmul_readvariableop_resource:2A
3sequential_1_output_biasadd_readvariableop_resource:
identityИв*sequential_1/conv2d/BiasAdd/ReadVariableOpв)sequential_1/conv2d/Conv2D/ReadVariableOpв,sequential_1/conv2d_1/BiasAdd/ReadVariableOpв+sequential_1/conv2d_1/Conv2D/ReadVariableOpв,sequential_1/conv2d_2/BiasAdd/ReadVariableOpв+sequential_1/conv2d_2/Conv2D/ReadVariableOpв,sequential_1/conv2d_3/BiasAdd/ReadVariableOpв+sequential_1/conv2d_3/Conv2D/ReadVariableOpв,sequential_1/conv2d_4/BiasAdd/ReadVariableOpв+sequential_1/conv2d_4/Conv2D/ReadVariableOpв)sequential_1/dense/BiasAdd/ReadVariableOpв(sequential_1/dense/MatMul/ReadVariableOpв+sequential_1/dense_1/BiasAdd/ReadVariableOpв*sequential_1/dense_1/MatMul/ReadVariableOpв+sequential_1/dense_2/BiasAdd/ReadVariableOpв*sequential_1/dense_2/MatMul/ReadVariableOpв+sequential_1/dense_3/BiasAdd/ReadVariableOpв*sequential_1/dense_3/MatMul/ReadVariableOpв*sequential_1/output/BiasAdd/ReadVariableOpв)sequential_1/output/MatMul/ReadVariableOpЗ
sequential_1/rescaling_1/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *БАА;2!
sequential_1/rescaling_1/Cast/xЛ
!sequential_1/rescaling_1/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!sequential_1/rescaling_1/Cast_1/x╣
sequential_1/rescaling_1/mulMulsequential_input(sequential_1/rescaling_1/Cast/x:output:0*
T0*/
_output_shapes
:         &2
sequential_1/rescaling_1/mul═
sequential_1/rescaling_1/addAddV2 sequential_1/rescaling_1/mul:z:0*sequential_1/rescaling_1/Cast_1/x:output:0*
T0*/
_output_shapes
:         &2
sequential_1/rescaling_1/add╤
)sequential_1/conv2d/Conv2D/ReadVariableOpReadVariableOp2sequential_1_conv2d_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02+
)sequential_1/conv2d/Conv2D/ReadVariableOp∙
sequential_1/conv2d/Conv2DConv2D sequential_1/rescaling_1/add:z:01sequential_1/conv2d/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &*
paddingSAME*
strides
2
sequential_1/conv2d/Conv2D╚
*sequential_1/conv2d/BiasAdd/ReadVariableOpReadVariableOp3sequential_1_conv2d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02,
*sequential_1/conv2d/BiasAdd/ReadVariableOp╪
sequential_1/conv2d/BiasAddBiasAdd#sequential_1/conv2d/Conv2D:output:02sequential_1/conv2d/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &2
sequential_1/conv2d/BiasAddЬ
sequential_1/conv2d/ReluRelu$sequential_1/conv2d/BiasAdd:output:0*
T0*/
_output_shapes
:         &2
sequential_1/conv2d/Reluш
"sequential_1/max_pooling2d/MaxPoolMaxPool&sequential_1/conv2d/Relu:activations:0*/
_output_shapes
:         *
ksize
*
paddingVALID*
strides
2$
"sequential_1/max_pooling2d/MaxPool╫
+sequential_1/conv2d_1/Conv2D/ReadVariableOpReadVariableOp4sequential_1_conv2d_1_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02-
+sequential_1/conv2d_1/Conv2D/ReadVariableOpК
sequential_1/conv2d_1/Conv2DConv2D+sequential_1/max_pooling2d/MaxPool:output:03sequential_1/conv2d_1/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         *
paddingSAME*
strides
2
sequential_1/conv2d_1/Conv2D╬
,sequential_1/conv2d_1/BiasAdd/ReadVariableOpReadVariableOp5sequential_1_conv2d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,sequential_1/conv2d_1/BiasAdd/ReadVariableOpр
sequential_1/conv2d_1/BiasAddBiasAdd%sequential_1/conv2d_1/Conv2D:output:04sequential_1/conv2d_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
sequential_1/conv2d_1/BiasAddв
sequential_1/conv2d_1/ReluRelu&sequential_1/conv2d_1/BiasAdd:output:0*
T0*/
_output_shapes
:         2
sequential_1/conv2d_1/Reluю
$sequential_1/max_pooling2d_1/MaxPoolMaxPool(sequential_1/conv2d_1/Relu:activations:0*/
_output_shapes
:         	*
ksize
*
paddingVALID*
strides
2&
$sequential_1/max_pooling2d_1/MaxPool╫
+sequential_1/conv2d_2/Conv2D/ReadVariableOpReadVariableOp4sequential_1_conv2d_2_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype02-
+sequential_1/conv2d_2/Conv2D/ReadVariableOpМ
sequential_1/conv2d_2/Conv2DConv2D-sequential_1/max_pooling2d_1/MaxPool:output:03sequential_1/conv2d_2/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 *
paddingSAME*
strides
2
sequential_1/conv2d_2/Conv2D╬
,sequential_1/conv2d_2/BiasAdd/ReadVariableOpReadVariableOp5sequential_1_conv2d_2_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02.
,sequential_1/conv2d_2/BiasAdd/ReadVariableOpр
sequential_1/conv2d_2/BiasAddBiasAdd%sequential_1/conv2d_2/Conv2D:output:04sequential_1/conv2d_2/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 2
sequential_1/conv2d_2/BiasAddв
sequential_1/conv2d_2/ReluRelu&sequential_1/conv2d_2/BiasAdd:output:0*
T0*/
_output_shapes
:         	 2
sequential_1/conv2d_2/Reluю
$sequential_1/max_pooling2d_2/MaxPoolMaxPool(sequential_1/conv2d_2/Relu:activations:0*/
_output_shapes
:          *
ksize
*
paddingVALID*
strides
2&
$sequential_1/max_pooling2d_2/MaxPool╫
+sequential_1/conv2d_3/Conv2D/ReadVariableOpReadVariableOp4sequential_1_conv2d_3_conv2d_readvariableop_resource*&
_output_shapes
: @*
dtype02-
+sequential_1/conv2d_3/Conv2D/ReadVariableOpМ
sequential_1/conv2d_3/Conv2DConv2D-sequential_1/max_pooling2d_2/MaxPool:output:03sequential_1/conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @*
paddingSAME*
strides
2
sequential_1/conv2d_3/Conv2D╬
,sequential_1/conv2d_3/BiasAdd/ReadVariableOpReadVariableOp5sequential_1_conv2d_3_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02.
,sequential_1/conv2d_3/BiasAdd/ReadVariableOpр
sequential_1/conv2d_3/BiasAddBiasAdd%sequential_1/conv2d_3/Conv2D:output:04sequential_1/conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @2
sequential_1/conv2d_3/BiasAddв
sequential_1/conv2d_3/ReluRelu&sequential_1/conv2d_3/BiasAdd:output:0*
T0*/
_output_shapes
:         @2
sequential_1/conv2d_3/Relu╫
+sequential_1/conv2d_4/Conv2D/ReadVariableOpReadVariableOp4sequential_1_conv2d_4_conv2d_readvariableop_resource*&
_output_shapes
:@P*
dtype02-
+sequential_1/conv2d_4/Conv2D/ReadVariableOpЗ
sequential_1/conv2d_4/Conv2DConv2D(sequential_1/conv2d_3/Relu:activations:03sequential_1/conv2d_4/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P*
paddingSAME*
strides
2
sequential_1/conv2d_4/Conv2D╬
,sequential_1/conv2d_4/BiasAdd/ReadVariableOpReadVariableOp5sequential_1_conv2d_4_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02.
,sequential_1/conv2d_4/BiasAdd/ReadVariableOpр
sequential_1/conv2d_4/BiasAddBiasAdd%sequential_1/conv2d_4/Conv2D:output:04sequential_1/conv2d_4/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P2
sequential_1/conv2d_4/BiasAddв
sequential_1/conv2d_4/ReluRelu&sequential_1/conv2d_4/BiasAdd:output:0*
T0*/
_output_shapes
:         P2
sequential_1/conv2d_4/Reluю
$sequential_1/max_pooling2d_3/MaxPoolMaxPool(sequential_1/conv2d_4/Relu:activations:0*/
_output_shapes
:         P*
ksize
*
paddingVALID*
strides
2&
$sequential_1/max_pooling2d_3/MaxPoolЙ
sequential_1/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"    а   2
sequential_1/flatten/Const╬
sequential_1/flatten/ReshapeReshape-sequential_1/max_pooling2d_3/MaxPool:output:0#sequential_1/flatten/Const:output:0*
T0*(
_output_shapes
:         а2
sequential_1/flatten/Reshape╚
(sequential_1/dense/MatMul/ReadVariableOpReadVariableOp1sequential_1_dense_matmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype02*
(sequential_1/dense/MatMul/ReadVariableOp╠
sequential_1/dense/MatMulMatMul%sequential_1/flatten/Reshape:output:00sequential_1/dense/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2
sequential_1/dense/MatMul╞
)sequential_1/dense/BiasAdd/ReadVariableOpReadVariableOp2sequential_1_dense_biasadd_readvariableop_resource*
_output_shapes	
:ю*
dtype02+
)sequential_1/dense/BiasAdd/ReadVariableOp╬
sequential_1/dense/BiasAddBiasAdd#sequential_1/dense/MatMul:product:01sequential_1/dense/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2
sequential_1/dense/BiasAddТ
sequential_1/dense/ReluRelu#sequential_1/dense/BiasAdd:output:0*
T0*(
_output_shapes
:         ю2
sequential_1/dense/Relu╬
*sequential_1/dense_1/MatMul/ReadVariableOpReadVariableOp3sequential_1_dense_1_matmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype02,
*sequential_1/dense_1/MatMul/ReadVariableOp╥
sequential_1/dense_1/MatMulMatMul%sequential_1/dense/Relu:activations:02sequential_1/dense_1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2
sequential_1/dense_1/MatMul╠
+sequential_1/dense_1/BiasAdd/ReadVariableOpReadVariableOp4sequential_1_dense_1_biasadd_readvariableop_resource*
_output_shapes	
:м*
dtype02-
+sequential_1/dense_1/BiasAdd/ReadVariableOp╓
sequential_1/dense_1/BiasAddBiasAdd%sequential_1/dense_1/MatMul:product:03sequential_1/dense_1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2
sequential_1/dense_1/BiasAddШ
sequential_1/dense_1/ReluRelu%sequential_1/dense_1/BiasAdd:output:0*
T0*(
_output_shapes
:         м2
sequential_1/dense_1/Reluж
sequential_1/dropout/IdentityIdentity'sequential_1/dense_1/Relu:activations:0*
T0*(
_output_shapes
:         м2
sequential_1/dropout/Identity═
*sequential_1/dense_2/MatMul/ReadVariableOpReadVariableOp3sequential_1_dense_2_matmul_readvariableop_resource*
_output_shapes
:	мd*
dtype02,
*sequential_1/dense_2/MatMul/ReadVariableOp╥
sequential_1/dense_2/MatMulMatMul&sequential_1/dropout/Identity:output:02sequential_1/dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2
sequential_1/dense_2/MatMul╦
+sequential_1/dense_2/BiasAdd/ReadVariableOpReadVariableOp4sequential_1_dense_2_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02-
+sequential_1/dense_2/BiasAdd/ReadVariableOp╒
sequential_1/dense_2/BiasAddBiasAdd%sequential_1/dense_2/MatMul:product:03sequential_1/dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2
sequential_1/dense_2/BiasAddЧ
sequential_1/dense_2/ReluRelu%sequential_1/dense_2/BiasAdd:output:0*
T0*'
_output_shapes
:         d2
sequential_1/dense_2/Relu╠
*sequential_1/dense_3/MatMul/ReadVariableOpReadVariableOp3sequential_1_dense_3_matmul_readvariableop_resource*
_output_shapes

:d2*
dtype02,
*sequential_1/dense_3/MatMul/ReadVariableOp╙
sequential_1/dense_3/MatMulMatMul'sequential_1/dense_2/Relu:activations:02sequential_1/dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22
sequential_1/dense_3/MatMul╦
+sequential_1/dense_3/BiasAdd/ReadVariableOpReadVariableOp4sequential_1_dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype02-
+sequential_1/dense_3/BiasAdd/ReadVariableOp╒
sequential_1/dense_3/BiasAddBiasAdd%sequential_1/dense_3/MatMul:product:03sequential_1/dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22
sequential_1/dense_3/BiasAddЧ
sequential_1/dense_3/ReluRelu%sequential_1/dense_3/BiasAdd:output:0*
T0*'
_output_shapes
:         22
sequential_1/dense_3/Reluй
sequential_1/dropout_1/IdentityIdentity'sequential_1/dense_3/Relu:activations:0*
T0*'
_output_shapes
:         22!
sequential_1/dropout_1/Identity╔
)sequential_1/output/MatMul/ReadVariableOpReadVariableOp2sequential_1_output_matmul_readvariableop_resource*
_output_shapes

:2*
dtype02+
)sequential_1/output/MatMul/ReadVariableOp╤
sequential_1/output/MatMulMatMul(sequential_1/dropout_1/Identity:output:01sequential_1/output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
sequential_1/output/MatMul╚
*sequential_1/output/BiasAdd/ReadVariableOpReadVariableOp3sequential_1_output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02,
*sequential_1/output/BiasAdd/ReadVariableOp╤
sequential_1/output/BiasAddBiasAdd$sequential_1/output/MatMul:product:02sequential_1/output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
sequential_1/output/BiasAddЭ
sequential_1/output/SoftmaxSoftmax$sequential_1/output/BiasAdd:output:0*
T0*'
_output_shapes
:         2
sequential_1/output/SoftmaxА
IdentityIdentity%sequential_1/output/Softmax:softmax:0^NoOp*
T0*'
_output_shapes
:         2

Identity▄
NoOpNoOp+^sequential_1/conv2d/BiasAdd/ReadVariableOp*^sequential_1/conv2d/Conv2D/ReadVariableOp-^sequential_1/conv2d_1/BiasAdd/ReadVariableOp,^sequential_1/conv2d_1/Conv2D/ReadVariableOp-^sequential_1/conv2d_2/BiasAdd/ReadVariableOp,^sequential_1/conv2d_2/Conv2D/ReadVariableOp-^sequential_1/conv2d_3/BiasAdd/ReadVariableOp,^sequential_1/conv2d_3/Conv2D/ReadVariableOp-^sequential_1/conv2d_4/BiasAdd/ReadVariableOp,^sequential_1/conv2d_4/Conv2D/ReadVariableOp*^sequential_1/dense/BiasAdd/ReadVariableOp)^sequential_1/dense/MatMul/ReadVariableOp,^sequential_1/dense_1/BiasAdd/ReadVariableOp+^sequential_1/dense_1/MatMul/ReadVariableOp,^sequential_1/dense_2/BiasAdd/ReadVariableOp+^sequential_1/dense_2/MatMul/ReadVariableOp,^sequential_1/dense_3/BiasAdd/ReadVariableOp+^sequential_1/dense_3/MatMul/ReadVariableOp+^sequential_1/output/BiasAdd/ReadVariableOp*^sequential_1/output/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*V
_input_shapesE
C:         &: : : : : : : : : : : : : : : : : : : : 2X
*sequential_1/conv2d/BiasAdd/ReadVariableOp*sequential_1/conv2d/BiasAdd/ReadVariableOp2V
)sequential_1/conv2d/Conv2D/ReadVariableOp)sequential_1/conv2d/Conv2D/ReadVariableOp2\
,sequential_1/conv2d_1/BiasAdd/ReadVariableOp,sequential_1/conv2d_1/BiasAdd/ReadVariableOp2Z
+sequential_1/conv2d_1/Conv2D/ReadVariableOp+sequential_1/conv2d_1/Conv2D/ReadVariableOp2\
,sequential_1/conv2d_2/BiasAdd/ReadVariableOp,sequential_1/conv2d_2/BiasAdd/ReadVariableOp2Z
+sequential_1/conv2d_2/Conv2D/ReadVariableOp+sequential_1/conv2d_2/Conv2D/ReadVariableOp2\
,sequential_1/conv2d_3/BiasAdd/ReadVariableOp,sequential_1/conv2d_3/BiasAdd/ReadVariableOp2Z
+sequential_1/conv2d_3/Conv2D/ReadVariableOp+sequential_1/conv2d_3/Conv2D/ReadVariableOp2\
,sequential_1/conv2d_4/BiasAdd/ReadVariableOp,sequential_1/conv2d_4/BiasAdd/ReadVariableOp2Z
+sequential_1/conv2d_4/Conv2D/ReadVariableOp+sequential_1/conv2d_4/Conv2D/ReadVariableOp2V
)sequential_1/dense/BiasAdd/ReadVariableOp)sequential_1/dense/BiasAdd/ReadVariableOp2T
(sequential_1/dense/MatMul/ReadVariableOp(sequential_1/dense/MatMul/ReadVariableOp2Z
+sequential_1/dense_1/BiasAdd/ReadVariableOp+sequential_1/dense_1/BiasAdd/ReadVariableOp2X
*sequential_1/dense_1/MatMul/ReadVariableOp*sequential_1/dense_1/MatMul/ReadVariableOp2Z
+sequential_1/dense_2/BiasAdd/ReadVariableOp+sequential_1/dense_2/BiasAdd/ReadVariableOp2X
*sequential_1/dense_2/MatMul/ReadVariableOp*sequential_1/dense_2/MatMul/ReadVariableOp2Z
+sequential_1/dense_3/BiasAdd/ReadVariableOp+sequential_1/dense_3/BiasAdd/ReadVariableOp2X
*sequential_1/dense_3/MatMul/ReadVariableOp*sequential_1/dense_3/MatMul/ReadVariableOp2X
*sequential_1/output/BiasAdd/ReadVariableOp*sequential_1/output/BiasAdd/ReadVariableOp2V
)sequential_1/output/MatMul/ReadVariableOp)sequential_1/output/MatMul/ReadVariableOp:a ]
/
_output_shapes
:         &
*
_user_specified_namesequential_input
╗
d
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_78284

inputs
identityТ
MaxPoolMaxPoolinputs*/
_output_shapes
:         *
ksize
*
paddingVALID*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:         2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ё
Х
'__inference_dense_2_layer_call_fn_78541

inputs
unknown:	мd
	unknown_0:d
identityИвStatefulPartitionedCallЄ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         d*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_765172
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         d2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         м: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
═
е
@__inference_dense_layer_call_and_return_conditional_losses_76464

inputs2
matmul_readvariableop_resource:
аю.
biasadd_readvariableop_resource:	ю
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpв.dense/kernel/Regularizer/Square/ReadVariableOpП
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:ю*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         ю2
Relu┴
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/muln
IdentityIdentityRelu:activations:0^NoOp*
T0*(
_output_shapes
:         ю2

Identity░
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp/^dense/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         а: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:         а
 
_user_specified_nameinputs
Ъ
b
F__inference_random_flip_layer_call_and_return_conditional_losses_78697

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
ф
^
B__inference_flatten_layer_call_and_return_conditional_losses_78435

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"    а   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:         а2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:         а2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         P:W S
/
_output_shapes
:         P
 
_user_specified_nameinputs
р
G
+__inference_random_flip_layer_call_fn_78686

inputs
identity╠
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_flip_layer_call_and_return_conditional_losses_758132
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
я
з
B__inference_dense_2_layer_call_and_return_conditional_losses_78558

inputs1
matmul_readvariableop_resource:	мd-
biasadd_readvariableop_resource:d
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpв0dense_2/kernel/Regularizer/Square/ReadVariableOpО
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	мd*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:         d2
Relu─
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mulm
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:         d2

Identity▓
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         м: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
╤
I
-__inference_max_pooling2d_layer_call_fn_78269

inputs
identityщ
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4                                    * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *Q
fLRJ
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_762332
PartitionedCallП
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
╒
K
/__inference_max_pooling2d_2_layer_call_fn_78349

inputs
identityы
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4                                    * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_762772
PartitionedCallП
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
│
м
*__inference_sequential_layer_call_fn_77941

inputs
unknown:	
	unknown_0:	
	unknown_1:	
identityИвStatefulPartitionedCallЕ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *N
fIRG
E__inference_sequential_layer_call_and_return_conditional_losses_761842
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         &2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*4
_input_shapes#
!:         &: : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
╣
a
E__inference_sequential_layer_call_and_return_conditional_losses_75828

inputs
identityф
random_flip/PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_flip_layer_call_and_return_conditional_losses_758132
random_flip/PartitionedCallО
random_rotation/PartitionedCallPartitionedCall$random_flip/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_random_rotation_layer_call_and_return_conditional_losses_758192!
random_rotation/PartitionedCallЖ
random_zoom/PartitionedCallPartitionedCall(random_rotation/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_zoom_layer_call_and_return_conditional_losses_758252
random_zoom/PartitionedCallА
IdentityIdentity$random_zoom/PartitionedCall:output:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
к
f
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_78359

inputs
identityн
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4                                    *
ksize
*
paddingVALID*
strides
2	
MaxPoolЗ
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
ї
м
__inference_loss_fn_0_78648K
7dense_kernel_regularizer_square_readvariableop_resource:
аю
identityИв.dense/kernel/Regularizer/Square/ReadVariableOp┌
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mulj
IdentityIdentity dense/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity
NoOpNoOp/^dense/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp
╝
E
)__inference_dropout_1_layer_call_fn_78595

inputs
identity┬
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *M
fHRF
D__inference_dropout_1_layer_call_and_return_conditional_losses_765512
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:         22

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:         2:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
ц
{
+__inference_random_flip_layer_call_fn_78693

inputs
unknown:	
identityИвStatefulPartitionedCallю
StatefulPartitionedCallStatefulPartitionedCallinputsunknown*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_flip_layer_call_and_return_conditional_losses_761592
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         &2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ю
b
F__inference_rescaling_1_layer_call_and_return_conditional_losses_76328

inputs
identityU
Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *БАА;2
Cast/xY
Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2

Cast_1/xd
mulMulinputsCast/x:output:0*
T0*/
_output_shapes
:         &2
muli
addAddV2mul:z:0Cast_1/x:output:0*
T0*/
_output_shapes
:         &2
addc
IdentityIdentityadd:z:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
є
`
B__inference_dropout_layer_call_and_return_conditional_losses_76498

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:         м2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:         м2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:         м:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
щ
№
C__inference_conv2d_1_layer_call_and_return_conditional_losses_76364

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         *
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:         
 
_user_specified_nameinputs
щ
№
C__inference_conv2d_2_layer_call_and_return_conditional_losses_76387

inputs8
conv2d_readvariableop_resource: -
biasadd_readvariableop_resource: 
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: *
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 *
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         	 2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         	 2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:         	
 
_user_specified_nameinputs
╩
C
'__inference_flatten_layer_call_fn_78429

inputs
identity┴
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         а* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_764452
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:         а2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         P:W S
/
_output_shapes
:         P
 
_user_specified_nameinputs
╜
f
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_76374

inputs
identityТ
MaxPoolMaxPoolinputs*/
_output_shapes
:         	*
ksize
*
paddingVALID*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:         	2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         :W S
/
_output_shapes
:         
 
_user_specified_nameinputs
·
┴
E__inference_sequential_layer_call_and_return_conditional_losses_76184

inputs
random_flip_76174:	#
random_rotation_76177:	
random_zoom_76180:	
identityИв#random_flip/StatefulPartitionedCallв'random_rotation/StatefulPartitionedCallв#random_zoom/StatefulPartitionedCallР
#random_flip/StatefulPartitionedCallStatefulPartitionedCallinputsrandom_flip_76174*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_flip_layer_call_and_return_conditional_losses_761592%
#random_flip/StatefulPartitionedCall╞
'random_rotation/StatefulPartitionedCallStatefulPartitionedCall,random_flip/StatefulPartitionedCall:output:0random_rotation_76177*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_random_rotation_layer_call_and_return_conditional_losses_760882)
'random_rotation/StatefulPartitionedCall║
#random_zoom/StatefulPartitionedCallStatefulPartitionedCall0random_rotation/StatefulPartitionedCall:output:0random_zoom_76180*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_zoom_layer_call_and_return_conditional_losses_759532%
#random_zoom/StatefulPartitionedCallП
IdentityIdentity,random_zoom/StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         &2

Identity─
NoOpNoOp$^random_flip/StatefulPartitionedCall(^random_rotation/StatefulPartitionedCall$^random_zoom/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*4
_input_shapes#
!:         &: : : 2J
#random_flip/StatefulPartitionedCall#random_flip/StatefulPartitionedCall2R
'random_rotation/StatefulPartitionedCall'random_rotation/StatefulPartitionedCall2J
#random_zoom/StatefulPartitionedCall#random_zoom/StatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
 
Q
*__inference_sequential_layer_call_fn_75831
random_flip_input
identity╓
PartitionedCallPartitionedCallrandom_flip_input*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *N
fIRG
E__inference_sequential_layer_call_and_return_conditional_losses_758282
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:b ^
/
_output_shapes
:         &
+
_user_specified_namerandom_flip_input
▐
F
*__inference_sequential_layer_call_fn_77930

inputs
identity╦
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *N
fIRG
E__inference_sequential_layer_call_and_return_conditional_losses_758282
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ш
Э
(__inference_conv2d_2_layer_call_fn_78333

inputs!
unknown: 
	unknown_0: 
identityИвStatefulPartitionedCall√
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_2_layer_call_and_return_conditional_losses_763872
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         	 2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         	
 
_user_specified_nameinputs
┌
l
E__inference_sequential_layer_call_and_return_conditional_losses_76211
random_flip_input
identityя
random_flip/PartitionedCallPartitionedCallrandom_flip_input*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_flip_layer_call_and_return_conditional_losses_758132
random_flip/PartitionedCallО
random_rotation/PartitionedCallPartitionedCall$random_flip/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_random_rotation_layer_call_and_return_conditional_losses_758192!
random_rotation/PartitionedCallЖ
random_zoom/PartitionedCallPartitionedCall(random_rotation/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_zoom_layer_call_and_return_conditional_losses_758252
random_zoom/PartitionedCallА
IdentityIdentity$random_zoom/PartitionedCall:output:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:b ^
/
_output_shapes
:         &
+
_user_specified_namerandom_flip_input
к
╕
#__inference_signature_wrapper_77313
sequential_input!
unknown:
	unknown_0:#
	unknown_1:
	unknown_2:#
	unknown_3: 
	unknown_4: #
	unknown_5: @
	unknown_6:@#
	unknown_7:@P
	unknown_8:P
	unknown_9:
аю

unknown_10:	ю

unknown_11:
юм

unknown_12:	м

unknown_13:	мd

unknown_14:d

unknown_15:d2

unknown_16:2

unknown_17:2

unknown_18:
identityИвStatefulPartitionedCall═
StatefulPartitionedCallStatefulPartitionedCallsequential_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8В *)
f$R"
 __inference__wrapped_model_758022
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*V
_input_shapesE
C:         &: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:a ]
/
_output_shapes
:         &
*
_user_specified_namesequential_input
р
G
+__inference_random_zoom_layer_call_fn_78898

inputs
identity╠
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_zoom_layer_call_and_return_conditional_losses_758252
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ш
Э
(__inference_conv2d_1_layer_call_fn_78293

inputs!
unknown:
	unknown_0:
identityИвStatefulPartitionedCall√
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_1_layer_call_and_return_conditional_losses_763642
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         
 
_user_specified_nameinputs
╜
f
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_78364

inputs
identityТ
MaxPoolMaxPoolinputs*/
_output_shapes
:          *
ksize
*
paddingVALID*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:          2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         	 :W S
/
_output_shapes
:         	 
 
_user_specified_nameinputs
Ш
Э
(__inference_conv2d_4_layer_call_fn_78393

inputs!
unknown:@P
	unknown_0:P
identityИвStatefulPartitionedCall√
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_764272
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         P2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         @: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         @
 
_user_specified_nameinputs
г
п
__inference_loss_fn_2_78670L
9dense_2_kernel_regularizer_square_readvariableop_resource:	мd
identityИв0dense_2/kernel/Regularizer/Square/ReadVariableOp▀
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_2_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mull
IdentityIdentity"dense_2/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

IdentityБ
NoOpNoOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp
к
c
D__inference_dropout_1_layer_call_and_return_conditional_losses_78617

inputs
identityИc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n█╢?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:         22
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape┤
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:         2*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЪЩЩ>2
dropout/GreaterEqual/y╛
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:         22
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:         22
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:         22
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:         22

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:         2:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
иа
├
J__inference_random_rotation_layer_call_and_return_conditional_losses_76088

inputs6
(stateful_uniform_rngreadandskip_resource:	
identityИвstateful_uniform/RngReadAndSkipD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliceБ
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:
¤        2
strided_slice_1/stackЕ
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ь
strided_slice_1StridedSliceShape:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1^
CastCaststrided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
CastБ
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice_2/stackЕ
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ь
strided_slice_2StridedSliceShape:output:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_2b
Cast_1Caststrided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
Cast_1~
stateful_uniform/shapePackstrided_slice:output:0*
N*
T0*
_output_shapes
:2
stateful_uniform/shapeq
stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *█╔└2
stateful_uniform/minq
stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *█╔@2
stateful_uniform/maxz
stateful_uniform/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
stateful_uniform/ConstЩ
stateful_uniform/ProdProdstateful_uniform/shape:output:0stateful_uniform/Const:output:0*
T0*
_output_shapes
: 2
stateful_uniform/Prodt
stateful_uniform/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2
stateful_uniform/Cast/xК
stateful_uniform/Cast_1Caststateful_uniform/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
stateful_uniform/Cast_1┘
stateful_uniform/RngReadAndSkipRngReadAndSkip(stateful_uniform_rngreadandskip_resource stateful_uniform/Cast/x:output:0stateful_uniform/Cast_1:y:0*
_output_shapes
:2!
stateful_uniform/RngReadAndSkipЦ
$stateful_uniform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2&
$stateful_uniform/strided_slice/stackЪ
&stateful_uniform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice/stack_1Ъ
&stateful_uniform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice/stack_2╬
stateful_uniform/strided_sliceStridedSlice'stateful_uniform/RngReadAndSkip:value:0-stateful_uniform/strided_slice/stack:output:0/stateful_uniform/strided_slice/stack_1:output:0/stateful_uniform/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2 
stateful_uniform/strided_sliceЩ
stateful_uniform/BitcastBitcast'stateful_uniform/strided_slice:output:0*
T0	*
_output_shapes
:*

type02
stateful_uniform/BitcastЪ
&stateful_uniform/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice_1/stackЮ
(stateful_uniform/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(stateful_uniform/strided_slice_1/stack_1Ю
(stateful_uniform/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(stateful_uniform/strided_slice_1/stack_2╞
 stateful_uniform/strided_slice_1StridedSlice'stateful_uniform/RngReadAndSkip:value:0/stateful_uniform/strided_slice_1/stack:output:01stateful_uniform/strided_slice_1/stack_1:output:01stateful_uniform/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2"
 stateful_uniform/strided_slice_1Я
stateful_uniform/Bitcast_1Bitcast)stateful_uniform/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02
stateful_uniform/Bitcast_1а
-stateful_uniform/StatelessRandomUniformV2/algConst*
_output_shapes
: *
dtype0*
value	B :2/
-stateful_uniform/StatelessRandomUniformV2/alg╕
)stateful_uniform/StatelessRandomUniformV2StatelessRandomUniformV2stateful_uniform/shape:output:0#stateful_uniform/Bitcast_1:output:0!stateful_uniform/Bitcast:output:06stateful_uniform/StatelessRandomUniformV2/alg:output:0*#
_output_shapes
:         2+
)stateful_uniform/StatelessRandomUniformV2Т
stateful_uniform/subSubstateful_uniform/max:output:0stateful_uniform/min:output:0*
T0*
_output_shapes
: 2
stateful_uniform/subп
stateful_uniform/mulMul2stateful_uniform/StatelessRandomUniformV2:output:0stateful_uniform/sub:z:0*
T0*#
_output_shapes
:         2
stateful_uniform/mulФ
stateful_uniformAddV2stateful_uniform/mul:z:0stateful_uniform/min:output:0*
T0*#
_output_shapes
:         2
stateful_uniforms
rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub/y~
rotation_matrix/subSub
Cast_1:y:0rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/subu
rotation_matrix/CosCosstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Cosw
rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_1/yД
rotation_matrix/sub_1Sub
Cast_1:y:0 rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_1У
rotation_matrix/mulMulrotation_matrix/Cos:y:0rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/mulu
rotation_matrix/SinSinstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Sinw
rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_2/yВ
rotation_matrix/sub_2SubCast:y:0 rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_2Ч
rotation_matrix/mul_1Mulrotation_matrix/Sin:y:0rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/mul_1Ч
rotation_matrix/sub_3Subrotation_matrix/mul:z:0rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/sub_3Ч
rotation_matrix/sub_4Subrotation_matrix/sub:z:0rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/sub_4{
rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
rotation_matrix/truediv/yк
rotation_matrix/truedivRealDivrotation_matrix/sub_4:z:0"rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:         2
rotation_matrix/truedivw
rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_5/yВ
rotation_matrix/sub_5SubCast:y:0 rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_5y
rotation_matrix/Sin_1Sinstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Sin_1w
rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_6/yД
rotation_matrix/sub_6Sub
Cast_1:y:0 rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_6Щ
rotation_matrix/mul_2Mulrotation_matrix/Sin_1:y:0rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/mul_2y
rotation_matrix/Cos_1Cosstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Cos_1w
rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_7/yВ
rotation_matrix/sub_7SubCast:y:0 rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_7Щ
rotation_matrix/mul_3Mulrotation_matrix/Cos_1:y:0rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/mul_3Ч
rotation_matrix/addAddV2rotation_matrix/mul_2:z:0rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/addЧ
rotation_matrix/sub_8Subrotation_matrix/sub_5:z:0rotation_matrix/add:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/sub_8
rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
rotation_matrix/truediv_1/y░
rotation_matrix/truediv_1RealDivrotation_matrix/sub_8:z:0$rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:         2
rotation_matrix/truediv_1r
rotation_matrix/ShapeShapestateful_uniform:z:0*
T0*
_output_shapes
:2
rotation_matrix/ShapeФ
#rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2%
#rotation_matrix/strided_slice/stackШ
%rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%rotation_matrix/strided_slice/stack_1Ш
%rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%rotation_matrix/strided_slice/stack_2┬
rotation_matrix/strided_sliceStridedSlicerotation_matrix/Shape:output:0,rotation_matrix/strided_slice/stack:output:0.rotation_matrix/strided_slice/stack_1:output:0.rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
rotation_matrix/strided_slicey
rotation_matrix/Cos_2Cosstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Cos_2Я
%rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_1/stackг
'rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_1/stack_1г
'rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_1/stack_2ў
rotation_matrix/strided_slice_1StridedSlicerotation_matrix/Cos_2:y:0.rotation_matrix/strided_slice_1/stack:output:00rotation_matrix/strided_slice_1/stack_1:output:00rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_1y
rotation_matrix/Sin_2Sinstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Sin_2Я
%rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_2/stackг
'rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_2/stack_1г
'rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_2/stack_2ў
rotation_matrix/strided_slice_2StridedSlicerotation_matrix/Sin_2:y:0.rotation_matrix/strided_slice_2/stack:output:00rotation_matrix/strided_slice_2/stack_1:output:00rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_2Н
rotation_matrix/NegNeg(rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:         2
rotation_matrix/NegЯ
%rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_3/stackг
'rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_3/stack_1г
'rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_3/stack_2∙
rotation_matrix/strided_slice_3StridedSlicerotation_matrix/truediv:z:0.rotation_matrix/strided_slice_3/stack:output:00rotation_matrix/strided_slice_3/stack_1:output:00rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_3y
rotation_matrix/Sin_3Sinstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Sin_3Я
%rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_4/stackг
'rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_4/stack_1г
'rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_4/stack_2ў
rotation_matrix/strided_slice_4StridedSlicerotation_matrix/Sin_3:y:0.rotation_matrix/strided_slice_4/stack:output:00rotation_matrix/strided_slice_4/stack_1:output:00rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_4y
rotation_matrix/Cos_3Cosstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Cos_3Я
%rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_5/stackг
'rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_5/stack_1г
'rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_5/stack_2ў
rotation_matrix/strided_slice_5StridedSlicerotation_matrix/Cos_3:y:0.rotation_matrix/strided_slice_5/stack:output:00rotation_matrix/strided_slice_5/stack_1:output:00rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_5Я
%rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_6/stackг
'rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_6/stack_1г
'rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_6/stack_2√
rotation_matrix/strided_slice_6StridedSlicerotation_matrix/truediv_1:z:0.rotation_matrix/strided_slice_6/stack:output:00rotation_matrix/strided_slice_6/stack_1:output:00rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_6|
rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
rotation_matrix/zeros/mul/yм
rotation_matrix/zeros/mulMul&rotation_matrix/strided_slice:output:0$rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/zeros/mul
rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
rotation_matrix/zeros/Less/yз
rotation_matrix/zeros/LessLessrotation_matrix/zeros/mul:z:0%rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/zeros/LessВ
rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2 
rotation_matrix/zeros/packed/1├
rotation_matrix/zeros/packedPack&rotation_matrix/strided_slice:output:0'rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
rotation_matrix/zeros/packed
rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rotation_matrix/zeros/Const╡
rotation_matrix/zerosFill%rotation_matrix/zeros/packed:output:0$rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:         2
rotation_matrix/zeros|
rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
rotation_matrix/concat/axisи
rotation_matrix/concatConcatV2(rotation_matrix/strided_slice_1:output:0rotation_matrix/Neg:y:0(rotation_matrix/strided_slice_3:output:0(rotation_matrix/strided_slice_4:output:0(rotation_matrix/strided_slice_5:output:0(rotation_matrix/strided_slice_6:output:0rotation_matrix/zeros:output:0$rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2
rotation_matrix/concatX
transform/ShapeShapeinputs*
T0*
_output_shapes
:2
transform/ShapeИ
transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2
transform/strided_slice/stackМ
transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2!
transform/strided_slice/stack_1М
transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2!
transform/strided_slice/stack_2К
transform/strided_sliceStridedSlicetransform/Shape:output:0&transform/strided_slice/stack:output:0(transform/strided_slice/stack_1:output:0(transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2
transform/strided_sliceq
transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2
transform/fill_value╟
$transform/ImageProjectiveTransformV3ImageProjectiveTransformV3inputsrotation_matrix/concat:output:0 transform/strided_slice:output:0transform/fill_value:output:0*/
_output_shapes
:         &*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2&
$transform/ImageProjectiveTransformV3Ь
IdentityIdentity9transform/ImageProjectiveTransformV3:transformed_images:0^NoOp*
T0*/
_output_shapes
:         &2

Identityp
NoOpNoOp ^stateful_uniform/RngReadAndSkip*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 2B
stateful_uniform/RngReadAndSkipstateful_uniform/RngReadAndSkip:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
и
d
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_76233

inputs
identityн
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4                                    *
ksize
*
paddingVALID*
strides
2	
MaxPoolЗ
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
Ю
f
J__inference_random_rotation_layer_call_and_return_conditional_losses_75819

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
щ
№
C__inference_conv2d_3_layer_call_and_return_conditional_losses_76410

inputs8
conv2d_readvariableop_resource: @-
biasadd_readvariableop_resource:@
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: @*
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @*
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         @2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         @2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:          : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:          
 
_user_specified_nameinputs
Ш
Э
(__inference_conv2d_3_layer_call_fn_78373

inputs!
unknown: @
	unknown_0:@
identityИвStatefulPartitionedCall√
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         @*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_764102
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         @2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:          : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:          
 
_user_specified_nameinputs
Ъ
b
F__inference_random_zoom_layer_call_and_return_conditional_losses_75825

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
╘
╖
*__inference_sequential_layer_call_fn_76204
random_flip_input
unknown:	
	unknown_0:	
	unknown_1:	
identityИвStatefulPartitionedCallР
StatefulPartitionedCallStatefulPartitionedCallrandom_flip_inputunknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *N
fIRG
E__inference_sequential_layer_call_and_return_conditional_losses_761842
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         &2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*4
_input_shapes#
!:         &: : : 22
StatefulPartitionedCallStatefulPartitionedCall:b ^
/
_output_shapes
:         &
+
_user_specified_namerandom_flip_input
ф
I
-__inference_max_pooling2d_layer_call_fn_78274

inputs
identity╬
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *Q
fLRJ
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_763512
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
щ
№
C__inference_conv2d_3_layer_call_and_return_conditional_losses_78384

inputs8
conv2d_readvariableop_resource: @-
biasadd_readvariableop_resource:@
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: @*
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @*
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         @2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         @2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:          : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:          
 
_user_specified_nameinputs
ё
b
D__inference_dropout_1_layer_call_and_return_conditional_losses_76551

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:         22

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:         22

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:         2:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
╜
f
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_76397

inputs
identityТ
MaxPoolMaxPoolinputs*/
_output_shapes
:          *
ksize
*
paddingVALID*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:          2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         	 :W S
/
_output_shapes
:         	 
 
_user_specified_nameinputs
ш
K
/__inference_max_pooling2d_1_layer_call_fn_78314

inputs
identity╨
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_763742
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         	2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         :W S
/
_output_shapes
:         
 
_user_specified_nameinputs
ш
K
/__inference_max_pooling2d_3_layer_call_fn_78414

inputs
identity╨
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_764372
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         P2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         P:W S
/
_output_shapes
:         P
 
_user_specified_nameinputs
═
е
@__inference_dense_layer_call_and_return_conditional_losses_78467

inputs2
matmul_readvariableop_resource:
аю.
biasadd_readvariableop_resource:	ю
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpв.dense/kernel/Regularizer/Square/ReadVariableOpП
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:ю*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         ю2
Relu┴
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/muln
IdentityIdentityRelu:activations:0^NoOp*
T0*(
_output_shapes
:         ю2

Identity░
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp/^dense/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         а: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:         а
 
_user_specified_nameinputs
∙
й
B__inference_dense_1_layer_call_and_return_conditional_losses_78499

inputs2
matmul_readvariableop_resource:
юм.
biasadd_readvariableop_resource:	м
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpв0dense_1/kernel/Regularizer/Square/ReadVariableOpП
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2
MatMulН
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:м*
dtype02
BiasAdd/ReadVariableOpВ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2	
BiasAddY
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:         м2
Relu┼
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/muln
IdentityIdentityRelu:activations:0^NoOp*
T0*(
_output_shapes
:         м2

Identity▓
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_1/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         ю: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:         ю
 
_user_specified_nameinputs
МЧ
╙
G__inference_sequential_1_layer_call_and_return_conditional_losses_77519

inputs?
%conv2d_conv2d_readvariableop_resource:4
&conv2d_biasadd_readvariableop_resource:A
'conv2d_1_conv2d_readvariableop_resource:6
(conv2d_1_biasadd_readvariableop_resource:A
'conv2d_2_conv2d_readvariableop_resource: 6
(conv2d_2_biasadd_readvariableop_resource: A
'conv2d_3_conv2d_readvariableop_resource: @6
(conv2d_3_biasadd_readvariableop_resource:@A
'conv2d_4_conv2d_readvariableop_resource:@P6
(conv2d_4_biasadd_readvariableop_resource:P8
$dense_matmul_readvariableop_resource:
аю4
%dense_biasadd_readvariableop_resource:	ю:
&dense_1_matmul_readvariableop_resource:
юм6
'dense_1_biasadd_readvariableop_resource:	м9
&dense_2_matmul_readvariableop_resource:	мd5
'dense_2_biasadd_readvariableop_resource:d8
&dense_3_matmul_readvariableop_resource:d25
'dense_3_biasadd_readvariableop_resource:27
%output_matmul_readvariableop_resource:24
&output_biasadd_readvariableop_resource:
identityИвconv2d/BiasAdd/ReadVariableOpвconv2d/Conv2D/ReadVariableOpвconv2d_1/BiasAdd/ReadVariableOpвconv2d_1/Conv2D/ReadVariableOpвconv2d_2/BiasAdd/ReadVariableOpвconv2d_2/Conv2D/ReadVariableOpвconv2d_3/BiasAdd/ReadVariableOpвconv2d_3/Conv2D/ReadVariableOpвconv2d_4/BiasAdd/ReadVariableOpвconv2d_4/Conv2D/ReadVariableOpвdense/BiasAdd/ReadVariableOpвdense/MatMul/ReadVariableOpв.dense/kernel/Regularizer/Square/ReadVariableOpвdense_1/BiasAdd/ReadVariableOpвdense_1/MatMul/ReadVariableOpв0dense_1/kernel/Regularizer/Square/ReadVariableOpвdense_2/BiasAdd/ReadVariableOpвdense_2/MatMul/ReadVariableOpв0dense_2/kernel/Regularizer/Square/ReadVariableOpвdense_3/BiasAdd/ReadVariableOpвdense_3/MatMul/ReadVariableOpв0dense_3/kernel/Regularizer/Square/ReadVariableOpвoutput/BiasAdd/ReadVariableOpвoutput/MatMul/ReadVariableOpm
rescaling_1/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *БАА;2
rescaling_1/Cast/xq
rescaling_1/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rescaling_1/Cast_1/xИ
rescaling_1/mulMulinputsrescaling_1/Cast/x:output:0*
T0*/
_output_shapes
:         &2
rescaling_1/mulЩ
rescaling_1/addAddV2rescaling_1/mul:z:0rescaling_1/Cast_1/x:output:0*
T0*/
_output_shapes
:         &2
rescaling_1/addк
conv2d/Conv2D/ReadVariableOpReadVariableOp%conv2d_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
conv2d/Conv2D/ReadVariableOp┼
conv2d/Conv2DConv2Drescaling_1/add:z:0$conv2d/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &*
paddingSAME*
strides
2
conv2d/Conv2Dб
conv2d/BiasAdd/ReadVariableOpReadVariableOp&conv2d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
conv2d/BiasAdd/ReadVariableOpд
conv2d/BiasAddBiasAddconv2d/Conv2D:output:0%conv2d/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &2
conv2d/BiasAddu
conv2d/ReluReluconv2d/BiasAdd:output:0*
T0*/
_output_shapes
:         &2
conv2d/Relu┴
max_pooling2d/MaxPoolMaxPoolconv2d/Relu:activations:0*/
_output_shapes
:         *
ksize
*
paddingVALID*
strides
2
max_pooling2d/MaxPool░
conv2d_1/Conv2D/ReadVariableOpReadVariableOp'conv2d_1_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02 
conv2d_1/Conv2D/ReadVariableOp╓
conv2d_1/Conv2DConv2Dmax_pooling2d/MaxPool:output:0&conv2d_1/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         *
paddingSAME*
strides
2
conv2d_1/Conv2Dз
conv2d_1/BiasAdd/ReadVariableOpReadVariableOp(conv2d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv2d_1/BiasAdd/ReadVariableOpм
conv2d_1/BiasAddBiasAddconv2d_1/Conv2D:output:0'conv2d_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
conv2d_1/BiasAdd{
conv2d_1/ReluReluconv2d_1/BiasAdd:output:0*
T0*/
_output_shapes
:         2
conv2d_1/Relu╟
max_pooling2d_1/MaxPoolMaxPoolconv2d_1/Relu:activations:0*/
_output_shapes
:         	*
ksize
*
paddingVALID*
strides
2
max_pooling2d_1/MaxPool░
conv2d_2/Conv2D/ReadVariableOpReadVariableOp'conv2d_2_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype02 
conv2d_2/Conv2D/ReadVariableOp╪
conv2d_2/Conv2DConv2D max_pooling2d_1/MaxPool:output:0&conv2d_2/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 *
paddingSAME*
strides
2
conv2d_2/Conv2Dз
conv2d_2/BiasAdd/ReadVariableOpReadVariableOp(conv2d_2_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02!
conv2d_2/BiasAdd/ReadVariableOpм
conv2d_2/BiasAddBiasAddconv2d_2/Conv2D:output:0'conv2d_2/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 2
conv2d_2/BiasAdd{
conv2d_2/ReluReluconv2d_2/BiasAdd:output:0*
T0*/
_output_shapes
:         	 2
conv2d_2/Relu╟
max_pooling2d_2/MaxPoolMaxPoolconv2d_2/Relu:activations:0*/
_output_shapes
:          *
ksize
*
paddingVALID*
strides
2
max_pooling2d_2/MaxPool░
conv2d_3/Conv2D/ReadVariableOpReadVariableOp'conv2d_3_conv2d_readvariableop_resource*&
_output_shapes
: @*
dtype02 
conv2d_3/Conv2D/ReadVariableOp╪
conv2d_3/Conv2DConv2D max_pooling2d_2/MaxPool:output:0&conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @*
paddingSAME*
strides
2
conv2d_3/Conv2Dз
conv2d_3/BiasAdd/ReadVariableOpReadVariableOp(conv2d_3_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02!
conv2d_3/BiasAdd/ReadVariableOpм
conv2d_3/BiasAddBiasAddconv2d_3/Conv2D:output:0'conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @2
conv2d_3/BiasAdd{
conv2d_3/ReluReluconv2d_3/BiasAdd:output:0*
T0*/
_output_shapes
:         @2
conv2d_3/Relu░
conv2d_4/Conv2D/ReadVariableOpReadVariableOp'conv2d_4_conv2d_readvariableop_resource*&
_output_shapes
:@P*
dtype02 
conv2d_4/Conv2D/ReadVariableOp╙
conv2d_4/Conv2DConv2Dconv2d_3/Relu:activations:0&conv2d_4/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P*
paddingSAME*
strides
2
conv2d_4/Conv2Dз
conv2d_4/BiasAdd/ReadVariableOpReadVariableOp(conv2d_4_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02!
conv2d_4/BiasAdd/ReadVariableOpм
conv2d_4/BiasAddBiasAddconv2d_4/Conv2D:output:0'conv2d_4/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P2
conv2d_4/BiasAdd{
conv2d_4/ReluReluconv2d_4/BiasAdd:output:0*
T0*/
_output_shapes
:         P2
conv2d_4/Relu╟
max_pooling2d_3/MaxPoolMaxPoolconv2d_4/Relu:activations:0*/
_output_shapes
:         P*
ksize
*
paddingVALID*
strides
2
max_pooling2d_3/MaxPoolo
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"    а   2
flatten/ConstЪ
flatten/ReshapeReshape max_pooling2d_3/MaxPool:output:0flatten/Const:output:0*
T0*(
_output_shapes
:         а2
flatten/Reshapeб
dense/MatMul/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype02
dense/MatMul/ReadVariableOpШ
dense/MatMulMatMulflatten/Reshape:output:0#dense/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2
dense/MatMulЯ
dense/BiasAdd/ReadVariableOpReadVariableOp%dense_biasadd_readvariableop_resource*
_output_shapes	
:ю*
dtype02
dense/BiasAdd/ReadVariableOpЪ
dense/BiasAddBiasAdddense/MatMul:product:0$dense/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2
dense/BiasAddk

dense/ReluReludense/BiasAdd:output:0*
T0*(
_output_shapes
:         ю2

dense/Reluз
dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype02
dense_1/MatMul/ReadVariableOpЮ
dense_1/MatMulMatMuldense/Relu:activations:0%dense_1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2
dense_1/MatMulе
dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes	
:м*
dtype02 
dense_1/BiasAdd/ReadVariableOpв
dense_1/BiasAddBiasAdddense_1/MatMul:product:0&dense_1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2
dense_1/BiasAddq
dense_1/ReluReludense_1/BiasAdd:output:0*
T0*(
_output_shapes
:         м2
dense_1/Relu
dropout/IdentityIdentitydense_1/Relu:activations:0*
T0*(
_output_shapes
:         м2
dropout/Identityж
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes
:	мd*
dtype02
dense_2/MatMul/ReadVariableOpЮ
dense_2/MatMulMatMuldropout/Identity:output:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2
dense_2/MatMulд
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02 
dense_2/BiasAdd/ReadVariableOpб
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2
dense_2/BiasAddp
dense_2/ReluReludense_2/BiasAdd:output:0*
T0*'
_output_shapes
:         d2
dense_2/Reluе
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

:d2*
dtype02
dense_3/MatMul/ReadVariableOpЯ
dense_3/MatMulMatMuldense_2/Relu:activations:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22
dense_3/MatMulд
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype02 
dense_3/BiasAdd/ReadVariableOpб
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22
dense_3/BiasAddp
dense_3/ReluReludense_3/BiasAdd:output:0*
T0*'
_output_shapes
:         22
dense_3/ReluВ
dropout_1/IdentityIdentitydense_3/Relu:activations:0*
T0*'
_output_shapes
:         22
dropout_1/Identityв
output/MatMul/ReadVariableOpReadVariableOp%output_matmul_readvariableop_resource*
_output_shapes

:2*
dtype02
output/MatMul/ReadVariableOpЭ
output/MatMulMatMuldropout_1/Identity:output:0$output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
output/MatMulб
output/BiasAdd/ReadVariableOpReadVariableOp&output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
output/BiasAdd/ReadVariableOpЭ
output/BiasAddBiasAddoutput/MatMul:product:0%output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
output/BiasAddv
output/SoftmaxSoftmaxoutput/BiasAdd:output:0*
T0*'
_output_shapes
:         2
output/Softmax╟
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul═
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul╠
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul╦
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/muls
IdentityIdentityoutput/Softmax:softmax:0^NoOp*
T0*'
_output_shapes
:         2

Identityв
NoOpNoOp^conv2d/BiasAdd/ReadVariableOp^conv2d/Conv2D/ReadVariableOp ^conv2d_1/BiasAdd/ReadVariableOp^conv2d_1/Conv2D/ReadVariableOp ^conv2d_2/BiasAdd/ReadVariableOp^conv2d_2/Conv2D/ReadVariableOp ^conv2d_3/BiasAdd/ReadVariableOp^conv2d_3/Conv2D/ReadVariableOp ^conv2d_4/BiasAdd/ReadVariableOp^conv2d_4/Conv2D/ReadVariableOp^dense/BiasAdd/ReadVariableOp^dense/MatMul/ReadVariableOp/^dense/kernel/Regularizer/Square/ReadVariableOp^dense_1/BiasAdd/ReadVariableOp^dense_1/MatMul/ReadVariableOp1^dense_1/kernel/Regularizer/Square/ReadVariableOp^dense_2/BiasAdd/ReadVariableOp^dense_2/MatMul/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp^output/BiasAdd/ReadVariableOp^output/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*V
_input_shapesE
C:         &: : : : : : : : : : : : : : : : : : : : 2>
conv2d/BiasAdd/ReadVariableOpconv2d/BiasAdd/ReadVariableOp2<
conv2d/Conv2D/ReadVariableOpconv2d/Conv2D/ReadVariableOp2B
conv2d_1/BiasAdd/ReadVariableOpconv2d_1/BiasAdd/ReadVariableOp2@
conv2d_1/Conv2D/ReadVariableOpconv2d_1/Conv2D/ReadVariableOp2B
conv2d_2/BiasAdd/ReadVariableOpconv2d_2/BiasAdd/ReadVariableOp2@
conv2d_2/Conv2D/ReadVariableOpconv2d_2/Conv2D/ReadVariableOp2B
conv2d_3/BiasAdd/ReadVariableOpconv2d_3/BiasAdd/ReadVariableOp2@
conv2d_3/Conv2D/ReadVariableOpconv2d_3/Conv2D/ReadVariableOp2B
conv2d_4/BiasAdd/ReadVariableOpconv2d_4/BiasAdd/ReadVariableOp2@
conv2d_4/Conv2D/ReadVariableOpconv2d_4/Conv2D/ReadVariableOp2<
dense/BiasAdd/ReadVariableOpdense/BiasAdd/ReadVariableOp2:
dense/MatMul/ReadVariableOpdense/MatMul/ReadVariableOp2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp2@
dense_1/BiasAdd/ReadVariableOpdense_1/BiasAdd/ReadVariableOp2>
dense_1/MatMul/ReadVariableOpdense_1/MatMul/ReadVariableOp2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp2@
dense_2/BiasAdd/ReadVariableOpdense_2/BiasAdd/ReadVariableOp2>
dense_2/MatMul/ReadVariableOpdense_2/MatMul/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2>
output/BiasAdd/ReadVariableOpoutput/BiasAdd/ReadVariableOp2<
output/MatMul/ReadVariableOpoutput/MatMul/ReadVariableOp:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
к
f
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_76299

inputs
identityн
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4                                    *
ksize
*
paddingVALID*
strides
2	
MaxPoolЗ
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
є
`
B__inference_dropout_layer_call_and_return_conditional_losses_78514

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:         м2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:         м2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:         м:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
╝
╖
,__inference_sequential_1_layer_call_fn_77358

inputs!
unknown:
	unknown_0:#
	unknown_1:
	unknown_2:#
	unknown_3: 
	unknown_4: #
	unknown_5: @
	unknown_6:@#
	unknown_7:@P
	unknown_8:P
	unknown_9:
аю

unknown_10:	ю

unknown_11:
юм

unknown_12:	м

unknown_13:	мd

unknown_14:d

unknown_15:d2

unknown_16:2

unknown_17:2

unknown_18:
identityИвStatefulPartitionedCallъ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8В *P
fKRI
G__inference_sequential_1_layer_call_and_return_conditional_losses_765952
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*V
_input_shapesE
C:         &: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
╨f
№
F__inference_random_flip_layer_call_and_return_conditional_losses_78755

inputs?
1stateful_uniform_full_int_rngreadandskip_resource:	
identityИв(stateful_uniform_full_int/RngReadAndSkipвOstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgвVstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterМ
stateful_uniform_full_int/shapeConst*
_output_shapes
:*
dtype0*
valueB:2!
stateful_uniform_full_int/shapeМ
stateful_uniform_full_int/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
stateful_uniform_full_int/Const╜
stateful_uniform_full_int/ProdProd(stateful_uniform_full_int/shape:output:0(stateful_uniform_full_int/Const:output:0*
T0*
_output_shapes
: 2 
stateful_uniform_full_int/ProdЖ
 stateful_uniform_full_int/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2"
 stateful_uniform_full_int/Cast/xе
 stateful_uniform_full_int/Cast_1Cast'stateful_uniform_full_int/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2"
 stateful_uniform_full_int/Cast_1Ж
(stateful_uniform_full_int/RngReadAndSkipRngReadAndSkip1stateful_uniform_full_int_rngreadandskip_resource)stateful_uniform_full_int/Cast/x:output:0$stateful_uniform_full_int/Cast_1:y:0*
_output_shapes
:2*
(stateful_uniform_full_int/RngReadAndSkipи
-stateful_uniform_full_int/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2/
-stateful_uniform_full_int/strided_slice/stackм
/stateful_uniform_full_int/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:21
/stateful_uniform_full_int/strided_slice/stack_1м
/stateful_uniform_full_int/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:21
/stateful_uniform_full_int/strided_slice/stack_2Д
'stateful_uniform_full_int/strided_sliceStridedSlice0stateful_uniform_full_int/RngReadAndSkip:value:06stateful_uniform_full_int/strided_slice/stack:output:08stateful_uniform_full_int/strided_slice/stack_1:output:08stateful_uniform_full_int/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2)
'stateful_uniform_full_int/strided_slice┤
!stateful_uniform_full_int/BitcastBitcast0stateful_uniform_full_int/strided_slice:output:0*
T0	*
_output_shapes
:*

type02#
!stateful_uniform_full_int/Bitcastм
/stateful_uniform_full_int/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:21
/stateful_uniform_full_int/strided_slice_1/stack░
1stateful_uniform_full_int/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:23
1stateful_uniform_full_int/strided_slice_1/stack_1░
1stateful_uniform_full_int/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:23
1stateful_uniform_full_int/strided_slice_1/stack_2№
)stateful_uniform_full_int/strided_slice_1StridedSlice0stateful_uniform_full_int/RngReadAndSkip:value:08stateful_uniform_full_int/strided_slice_1/stack:output:0:stateful_uniform_full_int/strided_slice_1/stack_1:output:0:stateful_uniform_full_int/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2+
)stateful_uniform_full_int/strided_slice_1║
#stateful_uniform_full_int/Bitcast_1Bitcast2stateful_uniform_full_int/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02%
#stateful_uniform_full_int/Bitcast_1А
stateful_uniform_full_int/algConst*
_output_shapes
: *
dtype0*
value	B :2
stateful_uniform_full_int/algо
stateful_uniform_full_intStatelessRandomUniformFullIntV2(stateful_uniform_full_int/shape:output:0,stateful_uniform_full_int/Bitcast_1:output:0*stateful_uniform_full_int/Bitcast:output:0&stateful_uniform_full_int/alg:output:0*
_output_shapes
:*
dtype0	2
stateful_uniform_full_intb

zeros_likeConst*
_output_shapes
:*
dtype0	*
valueB	R 2

zeros_likeБ
stackPack"stateful_uniform_full_int:output:0zeros_like:output:0*
N*
T0	*
_output_shapes

:2
stack{
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        2
strided_slice/stack
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       2
strided_slice/stack_1
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2
strided_slice/stack_2И
strided_sliceStridedSlicestack:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask*
end_mask*
shrink_axis_mask2
strided_slice╙
3stateless_random_flip_left_right/control_dependencyIdentityinputs*
T0*
_class
loc:@inputs*/
_output_shapes
:         &25
3stateless_random_flip_left_right/control_dependency╝
&stateless_random_flip_left_right/ShapeShape<stateless_random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:2(
&stateless_random_flip_left_right/Shape╢
4stateless_random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 26
4stateless_random_flip_left_right/strided_slice/stack║
6stateless_random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:28
6stateless_random_flip_left_right/strided_slice/stack_1║
6stateless_random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:28
6stateless_random_flip_left_right/strided_slice/stack_2и
.stateless_random_flip_left_right/strided_sliceStridedSlice/stateless_random_flip_left_right/Shape:output:0=stateless_random_flip_left_right/strided_slice/stack:output:0?stateless_random_flip_left_right/strided_slice/stack_1:output:0?stateless_random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask20
.stateless_random_flip_left_right/strided_sliceё
?stateless_random_flip_left_right/stateless_random_uniform/shapePack7stateless_random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2A
?stateless_random_flip_left_right/stateless_random_uniform/shape├
=stateless_random_flip_left_right/stateless_random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2?
=stateless_random_flip_left_right/stateless_random_uniform/min├
=stateless_random_flip_left_right/stateless_random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2?
=stateless_random_flip_left_right/stateless_random_uniform/maxК
Vstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterStatelessRandomGetKeyCounterstrided_slice:output:0* 
_output_shapes
::2X
Vstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterм
Ostateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgStatelessRandomGetAlgW^stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter*
_output_shapes
: 2Q
Ostateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg╩
Rstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2StatelessRandomUniformV2Hstateless_random_flip_left_right/stateless_random_uniform/shape:output:0\stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:key:0`stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:counter:0Ustateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg:alg:0*#
_output_shapes
:         2T
Rstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2╢
=stateless_random_flip_left_right/stateless_random_uniform/subSubFstateless_random_flip_left_right/stateless_random_uniform/max:output:0Fstateless_random_flip_left_right/stateless_random_uniform/min:output:0*
T0*
_output_shapes
: 2?
=stateless_random_flip_left_right/stateless_random_uniform/sub╙
=stateless_random_flip_left_right/stateless_random_uniform/mulMul[stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2:output:0Astateless_random_flip_left_right/stateless_random_uniform/sub:z:0*
T0*#
_output_shapes
:         2?
=stateless_random_flip_left_right/stateless_random_uniform/mul╕
9stateless_random_flip_left_right/stateless_random_uniformAddV2Astateless_random_flip_left_right/stateless_random_uniform/mul:z:0Fstateless_random_flip_left_right/stateless_random_uniform/min:output:0*
T0*#
_output_shapes
:         2;
9stateless_random_flip_left_right/stateless_random_uniformж
0stateless_random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :22
0stateless_random_flip_left_right/Reshape/shape/1ж
0stateless_random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :22
0stateless_random_flip_left_right/Reshape/shape/2ж
0stateless_random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :22
0stateless_random_flip_left_right/Reshape/shape/3А
.stateless_random_flip_left_right/Reshape/shapePack7stateless_random_flip_left_right/strided_slice:output:09stateless_random_flip_left_right/Reshape/shape/1:output:09stateless_random_flip_left_right/Reshape/shape/2:output:09stateless_random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:20
.stateless_random_flip_left_right/Reshape/shapeС
(stateless_random_flip_left_right/ReshapeReshape=stateless_random_flip_left_right/stateless_random_uniform:z:07stateless_random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:         2*
(stateless_random_flip_left_right/Reshape╞
&stateless_random_flip_left_right/RoundRound1stateless_random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:         2(
&stateless_random_flip_left_right/Roundм
/stateless_random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:21
/stateless_random_flip_left_right/ReverseV2/axisЧ
*stateless_random_flip_left_right/ReverseV2	ReverseV2<stateless_random_flip_left_right/control_dependency:output:08stateless_random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:         &2,
*stateless_random_flip_left_right/ReverseV2ю
$stateless_random_flip_left_right/mulMul*stateless_random_flip_left_right/Round:y:03stateless_random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:         &2&
$stateless_random_flip_left_right/mulХ
&stateless_random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2(
&stateless_random_flip_left_right/sub/xъ
$stateless_random_flip_left_right/subSub/stateless_random_flip_left_right/sub/x:output:0*stateless_random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:         2&
$stateless_random_flip_left_right/sub∙
&stateless_random_flip_left_right/mul_1Mul(stateless_random_flip_left_right/sub:z:0<stateless_random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:         &2(
&stateless_random_flip_left_right/mul_1х
$stateless_random_flip_left_right/addAddV2(stateless_random_flip_left_right/mul:z:0*stateless_random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:         &2&
$stateless_random_flip_left_right/addЛ
IdentityIdentity(stateless_random_flip_left_right/add:z:0^NoOp*
T0*/
_output_shapes
:         &2

Identityд
NoOpNoOp)^stateful_uniform_full_int/RngReadAndSkipP^stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgW^stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 2T
(stateful_uniform_full_int/RngReadAndSkip(stateful_uniform_full_int/RngReadAndSkip2в
Ostateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgOstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg2░
Vstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterVstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
а
о
__inference_loss_fn_3_78681K
9dense_3_kernel_regularizer_square_readvariableop_resource:d2
identityИв0dense_3/kernel/Regularizer/Square/ReadVariableOp▐
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_3_kernel_regularizer_square_readvariableop_resource*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mull
IdentityIdentity"dense_3/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

IdentityБ
NoOpNoOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp
┌Є
Ч
E__inference_sequential_layer_call_and_return_conditional_losses_78231

inputsK
=random_flip_stateful_uniform_full_int_rngreadandskip_resource:	F
8random_rotation_stateful_uniform_rngreadandskip_resource:	B
4random_zoom_stateful_uniform_rngreadandskip_resource:	
identityИв4random_flip/stateful_uniform_full_int/RngReadAndSkipв[random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgвbrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterв/random_rotation/stateful_uniform/RngReadAndSkipв+random_zoom/stateful_uniform/RngReadAndSkipд
+random_flip/stateful_uniform_full_int/shapeConst*
_output_shapes
:*
dtype0*
valueB:2-
+random_flip/stateful_uniform_full_int/shapeд
+random_flip/stateful_uniform_full_int/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2-
+random_flip/stateful_uniform_full_int/Constэ
*random_flip/stateful_uniform_full_int/ProdProd4random_flip/stateful_uniform_full_int/shape:output:04random_flip/stateful_uniform_full_int/Const:output:0*
T0*
_output_shapes
: 2,
*random_flip/stateful_uniform_full_int/ProdЮ
,random_flip/stateful_uniform_full_int/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2.
,random_flip/stateful_uniform_full_int/Cast/x╔
,random_flip/stateful_uniform_full_int/Cast_1Cast3random_flip/stateful_uniform_full_int/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2.
,random_flip/stateful_uniform_full_int/Cast_1┬
4random_flip/stateful_uniform_full_int/RngReadAndSkipRngReadAndSkip=random_flip_stateful_uniform_full_int_rngreadandskip_resource5random_flip/stateful_uniform_full_int/Cast/x:output:00random_flip/stateful_uniform_full_int/Cast_1:y:0*
_output_shapes
:26
4random_flip/stateful_uniform_full_int/RngReadAndSkip└
9random_flip/stateful_uniform_full_int/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2;
9random_flip/stateful_uniform_full_int/strided_slice/stack─
;random_flip/stateful_uniform_full_int/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2=
;random_flip/stateful_uniform_full_int/strided_slice/stack_1─
;random_flip/stateful_uniform_full_int/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2=
;random_flip/stateful_uniform_full_int/strided_slice/stack_2╠
3random_flip/stateful_uniform_full_int/strided_sliceStridedSlice<random_flip/stateful_uniform_full_int/RngReadAndSkip:value:0Brandom_flip/stateful_uniform_full_int/strided_slice/stack:output:0Drandom_flip/stateful_uniform_full_int/strided_slice/stack_1:output:0Drandom_flip/stateful_uniform_full_int/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask25
3random_flip/stateful_uniform_full_int/strided_slice╪
-random_flip/stateful_uniform_full_int/BitcastBitcast<random_flip/stateful_uniform_full_int/strided_slice:output:0*
T0	*
_output_shapes
:*

type02/
-random_flip/stateful_uniform_full_int/Bitcast─
;random_flip/stateful_uniform_full_int/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2=
;random_flip/stateful_uniform_full_int/strided_slice_1/stack╚
=random_flip/stateful_uniform_full_int/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2?
=random_flip/stateful_uniform_full_int/strided_slice_1/stack_1╚
=random_flip/stateful_uniform_full_int/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2?
=random_flip/stateful_uniform_full_int/strided_slice_1/stack_2─
5random_flip/stateful_uniform_full_int/strided_slice_1StridedSlice<random_flip/stateful_uniform_full_int/RngReadAndSkip:value:0Drandom_flip/stateful_uniform_full_int/strided_slice_1/stack:output:0Frandom_flip/stateful_uniform_full_int/strided_slice_1/stack_1:output:0Frandom_flip/stateful_uniform_full_int/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:27
5random_flip/stateful_uniform_full_int/strided_slice_1▐
/random_flip/stateful_uniform_full_int/Bitcast_1Bitcast>random_flip/stateful_uniform_full_int/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type021
/random_flip/stateful_uniform_full_int/Bitcast_1Ш
)random_flip/stateful_uniform_full_int/algConst*
_output_shapes
: *
dtype0*
value	B :2+
)random_flip/stateful_uniform_full_int/algЎ
%random_flip/stateful_uniform_full_intStatelessRandomUniformFullIntV24random_flip/stateful_uniform_full_int/shape:output:08random_flip/stateful_uniform_full_int/Bitcast_1:output:06random_flip/stateful_uniform_full_int/Bitcast:output:02random_flip/stateful_uniform_full_int/alg:output:0*
_output_shapes
:*
dtype0	2'
%random_flip/stateful_uniform_full_intz
random_flip/zeros_likeConst*
_output_shapes
:*
dtype0	*
valueB	R 2
random_flip/zeros_like▒
random_flip/stackPack.random_flip/stateful_uniform_full_int:output:0random_flip/zeros_like:output:0*
N*
T0	*
_output_shapes

:2
random_flip/stackУ
random_flip/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        2!
random_flip/strided_slice/stackЧ
!random_flip/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       2#
!random_flip/strided_slice/stack_1Ч
!random_flip/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2#
!random_flip/strided_slice/stack_2╨
random_flip/strided_sliceStridedSlicerandom_flip/stack:output:0(random_flip/strided_slice/stack:output:0*random_flip/strided_slice/stack_1:output:0*random_flip/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask*
end_mask*
shrink_axis_mask2
random_flip/strided_sliceы
?random_flip/stateless_random_flip_left_right/control_dependencyIdentityinputs*
T0*
_class
loc:@inputs*/
_output_shapes
:         &2A
?random_flip/stateless_random_flip_left_right/control_dependencyр
2random_flip/stateless_random_flip_left_right/ShapeShapeHrandom_flip/stateless_random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:24
2random_flip/stateless_random_flip_left_right/Shape╬
@random_flip/stateless_random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2B
@random_flip/stateless_random_flip_left_right/strided_slice/stack╥
Brandom_flip/stateless_random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2D
Brandom_flip/stateless_random_flip_left_right/strided_slice/stack_1╥
Brandom_flip/stateless_random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2D
Brandom_flip/stateless_random_flip_left_right/strided_slice/stack_2Ё
:random_flip/stateless_random_flip_left_right/strided_sliceStridedSlice;random_flip/stateless_random_flip_left_right/Shape:output:0Irandom_flip/stateless_random_flip_left_right/strided_slice/stack:output:0Krandom_flip/stateless_random_flip_left_right/strided_slice/stack_1:output:0Krandom_flip/stateless_random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2<
:random_flip/stateless_random_flip_left_right/strided_sliceХ
Krandom_flip/stateless_random_flip_left_right/stateless_random_uniform/shapePackCrandom_flip/stateless_random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2M
Krandom_flip/stateless_random_flip_left_right/stateless_random_uniform/shape█
Irandom_flip/stateless_random_flip_left_right/stateless_random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2K
Irandom_flip/stateless_random_flip_left_right/stateless_random_uniform/min█
Irandom_flip/stateless_random_flip_left_right/stateless_random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2K
Irandom_flip/stateless_random_flip_left_right/stateless_random_uniform/maxо
brandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterStatelessRandomGetKeyCounter"random_flip/strided_slice:output:0* 
_output_shapes
::2d
brandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter╨
[random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgStatelessRandomGetAlgc^random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter*
_output_shapes
: 2]
[random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgТ
^random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2StatelessRandomUniformV2Trandom_flip/stateless_random_flip_left_right/stateless_random_uniform/shape:output:0hrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:key:0lrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:counter:0arandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg:alg:0*#
_output_shapes
:         2`
^random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2ц
Irandom_flip/stateless_random_flip_left_right/stateless_random_uniform/subSubRrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/max:output:0Rrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/min:output:0*
T0*
_output_shapes
: 2K
Irandom_flip/stateless_random_flip_left_right/stateless_random_uniform/subГ
Irandom_flip/stateless_random_flip_left_right/stateless_random_uniform/mulMulgrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2:output:0Mrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/sub:z:0*
T0*#
_output_shapes
:         2K
Irandom_flip/stateless_random_flip_left_right/stateless_random_uniform/mulш
Erandom_flip/stateless_random_flip_left_right/stateless_random_uniformAddV2Mrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/mul:z:0Rrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/min:output:0*
T0*#
_output_shapes
:         2G
Erandom_flip/stateless_random_flip_left_right/stateless_random_uniform╛
<random_flip/stateless_random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2>
<random_flip/stateless_random_flip_left_right/Reshape/shape/1╛
<random_flip/stateless_random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2>
<random_flip/stateless_random_flip_left_right/Reshape/shape/2╛
<random_flip/stateless_random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2>
<random_flip/stateless_random_flip_left_right/Reshape/shape/3╚
:random_flip/stateless_random_flip_left_right/Reshape/shapePackCrandom_flip/stateless_random_flip_left_right/strided_slice:output:0Erandom_flip/stateless_random_flip_left_right/Reshape/shape/1:output:0Erandom_flip/stateless_random_flip_left_right/Reshape/shape/2:output:0Erandom_flip/stateless_random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2<
:random_flip/stateless_random_flip_left_right/Reshape/shape┴
4random_flip/stateless_random_flip_left_right/ReshapeReshapeIrandom_flip/stateless_random_flip_left_right/stateless_random_uniform:z:0Crandom_flip/stateless_random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:         26
4random_flip/stateless_random_flip_left_right/Reshapeъ
2random_flip/stateless_random_flip_left_right/RoundRound=random_flip/stateless_random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:         24
2random_flip/stateless_random_flip_left_right/Round─
;random_flip/stateless_random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:2=
;random_flip/stateless_random_flip_left_right/ReverseV2/axis╟
6random_flip/stateless_random_flip_left_right/ReverseV2	ReverseV2Hrandom_flip/stateless_random_flip_left_right/control_dependency:output:0Drandom_flip/stateless_random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:         &28
6random_flip/stateless_random_flip_left_right/ReverseV2Ю
0random_flip/stateless_random_flip_left_right/mulMul6random_flip/stateless_random_flip_left_right/Round:y:0?random_flip/stateless_random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:         &22
0random_flip/stateless_random_flip_left_right/mulн
2random_flip/stateless_random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?24
2random_flip/stateless_random_flip_left_right/sub/xЪ
0random_flip/stateless_random_flip_left_right/subSub;random_flip/stateless_random_flip_left_right/sub/x:output:06random_flip/stateless_random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:         22
0random_flip/stateless_random_flip_left_right/subй
2random_flip/stateless_random_flip_left_right/mul_1Mul4random_flip/stateless_random_flip_left_right/sub:z:0Hrandom_flip/stateless_random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:         &24
2random_flip/stateless_random_flip_left_right/mul_1Х
0random_flip/stateless_random_flip_left_right/addAddV24random_flip/stateless_random_flip_left_right/mul:z:06random_flip/stateless_random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:         &22
0random_flip/stateless_random_flip_left_right/addТ
random_rotation/ShapeShape4random_flip/stateless_random_flip_left_right/add:z:0*
T0*
_output_shapes
:2
random_rotation/ShapeФ
#random_rotation/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2%
#random_rotation/strided_slice/stackШ
%random_rotation/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_rotation/strided_slice/stack_1Ш
%random_rotation/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%random_rotation/strided_slice/stack_2┬
random_rotation/strided_sliceStridedSlicerandom_rotation/Shape:output:0,random_rotation/strided_slice/stack:output:0.random_rotation/strided_slice/stack_1:output:0.random_rotation/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_rotation/strided_sliceб
%random_rotation/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:
¤        2'
%random_rotation/strided_slice_1/stackе
'random_rotation/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
■        2)
'random_rotation/strided_slice_1/stack_1Ь
'random_rotation/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation/strided_slice_1/stack_2╠
random_rotation/strided_slice_1StridedSlicerandom_rotation/Shape:output:0.random_rotation/strided_slice_1/stack:output:00random_rotation/strided_slice_1/stack_1:output:00random_rotation/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2!
random_rotation/strided_slice_1О
random_rotation/CastCast(random_rotation/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_rotation/Castб
%random_rotation/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2'
%random_rotation/strided_slice_2/stackе
'random_rotation/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2)
'random_rotation/strided_slice_2/stack_1Ь
'random_rotation/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2)
'random_rotation/strided_slice_2/stack_2╠
random_rotation/strided_slice_2StridedSlicerandom_rotation/Shape:output:0.random_rotation/strided_slice_2/stack:output:00random_rotation/strided_slice_2/stack_1:output:00random_rotation/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2!
random_rotation/strided_slice_2Т
random_rotation/Cast_1Cast(random_rotation/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_rotation/Cast_1о
&random_rotation/stateful_uniform/shapePack&random_rotation/strided_slice:output:0*
N*
T0*
_output_shapes
:2(
&random_rotation/stateful_uniform/shapeС
$random_rotation/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *█╔└2&
$random_rotation/stateful_uniform/minС
$random_rotation/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *█╔@2&
$random_rotation/stateful_uniform/maxЪ
&random_rotation/stateful_uniform/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2(
&random_rotation/stateful_uniform/Const┘
%random_rotation/stateful_uniform/ProdProd/random_rotation/stateful_uniform/shape:output:0/random_rotation/stateful_uniform/Const:output:0*
T0*
_output_shapes
: 2'
%random_rotation/stateful_uniform/ProdФ
'random_rotation/stateful_uniform/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2)
'random_rotation/stateful_uniform/Cast/x║
'random_rotation/stateful_uniform/Cast_1Cast.random_rotation/stateful_uniform/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2)
'random_rotation/stateful_uniform/Cast_1й
/random_rotation/stateful_uniform/RngReadAndSkipRngReadAndSkip8random_rotation_stateful_uniform_rngreadandskip_resource0random_rotation/stateful_uniform/Cast/x:output:0+random_rotation/stateful_uniform/Cast_1:y:0*
_output_shapes
:21
/random_rotation/stateful_uniform/RngReadAndSkip╢
4random_rotation/stateful_uniform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 26
4random_rotation/stateful_uniform/strided_slice/stack║
6random_rotation/stateful_uniform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:28
6random_rotation/stateful_uniform/strided_slice/stack_1║
6random_rotation/stateful_uniform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:28
6random_rotation/stateful_uniform/strided_slice/stack_2о
.random_rotation/stateful_uniform/strided_sliceStridedSlice7random_rotation/stateful_uniform/RngReadAndSkip:value:0=random_rotation/stateful_uniform/strided_slice/stack:output:0?random_rotation/stateful_uniform/strided_slice/stack_1:output:0?random_rotation/stateful_uniform/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask20
.random_rotation/stateful_uniform/strided_slice╔
(random_rotation/stateful_uniform/BitcastBitcast7random_rotation/stateful_uniform/strided_slice:output:0*
T0	*
_output_shapes
:*

type02*
(random_rotation/stateful_uniform/Bitcast║
6random_rotation/stateful_uniform/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:28
6random_rotation/stateful_uniform/strided_slice_1/stack╛
8random_rotation/stateful_uniform/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2:
8random_rotation/stateful_uniform/strided_slice_1/stack_1╛
8random_rotation/stateful_uniform/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2:
8random_rotation/stateful_uniform/strided_slice_1/stack_2ж
0random_rotation/stateful_uniform/strided_slice_1StridedSlice7random_rotation/stateful_uniform/RngReadAndSkip:value:0?random_rotation/stateful_uniform/strided_slice_1/stack:output:0Arandom_rotation/stateful_uniform/strided_slice_1/stack_1:output:0Arandom_rotation/stateful_uniform/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:22
0random_rotation/stateful_uniform/strided_slice_1╧
*random_rotation/stateful_uniform/Bitcast_1Bitcast9random_rotation/stateful_uniform/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02,
*random_rotation/stateful_uniform/Bitcast_1└
=random_rotation/stateful_uniform/StatelessRandomUniformV2/algConst*
_output_shapes
: *
dtype0*
value	B :2?
=random_rotation/stateful_uniform/StatelessRandomUniformV2/algШ
9random_rotation/stateful_uniform/StatelessRandomUniformV2StatelessRandomUniformV2/random_rotation/stateful_uniform/shape:output:03random_rotation/stateful_uniform/Bitcast_1:output:01random_rotation/stateful_uniform/Bitcast:output:0Frandom_rotation/stateful_uniform/StatelessRandomUniformV2/alg:output:0*#
_output_shapes
:         2;
9random_rotation/stateful_uniform/StatelessRandomUniformV2╥
$random_rotation/stateful_uniform/subSub-random_rotation/stateful_uniform/max:output:0-random_rotation/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2&
$random_rotation/stateful_uniform/subя
$random_rotation/stateful_uniform/mulMulBrandom_rotation/stateful_uniform/StatelessRandomUniformV2:output:0(random_rotation/stateful_uniform/sub:z:0*
T0*#
_output_shapes
:         2&
$random_rotation/stateful_uniform/mul╘
 random_rotation/stateful_uniformAddV2(random_rotation/stateful_uniform/mul:z:0-random_rotation/stateful_uniform/min:output:0*
T0*#
_output_shapes
:         2"
 random_rotation/stateful_uniformУ
%random_rotation/rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2'
%random_rotation/rotation_matrix/sub/y╛
#random_rotation/rotation_matrix/subSubrandom_rotation/Cast_1:y:0.random_rotation/rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 2%
#random_rotation/rotation_matrix/subе
#random_rotation/rotation_matrix/CosCos$random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         2%
#random_rotation/rotation_matrix/CosЧ
'random_rotation/rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2)
'random_rotation/rotation_matrix/sub_1/y─
%random_rotation/rotation_matrix/sub_1Subrandom_rotation/Cast_1:y:00random_rotation/rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 2'
%random_rotation/rotation_matrix/sub_1╙
#random_rotation/rotation_matrix/mulMul'random_rotation/rotation_matrix/Cos:y:0)random_rotation/rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:         2%
#random_rotation/rotation_matrix/mulе
#random_rotation/rotation_matrix/SinSin$random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         2%
#random_rotation/rotation_matrix/SinЧ
'random_rotation/rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2)
'random_rotation/rotation_matrix/sub_2/y┬
%random_rotation/rotation_matrix/sub_2Subrandom_rotation/Cast:y:00random_rotation/rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2'
%random_rotation/rotation_matrix/sub_2╫
%random_rotation/rotation_matrix/mul_1Mul'random_rotation/rotation_matrix/Sin:y:0)random_rotation/rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/mul_1╫
%random_rotation/rotation_matrix/sub_3Sub'random_rotation/rotation_matrix/mul:z:0)random_rotation/rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/sub_3╫
%random_rotation/rotation_matrix/sub_4Sub'random_rotation/rotation_matrix/sub:z:0)random_rotation/rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/sub_4Ы
)random_rotation/rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2+
)random_rotation/rotation_matrix/truediv/yъ
'random_rotation/rotation_matrix/truedivRealDiv)random_rotation/rotation_matrix/sub_4:z:02random_rotation/rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:         2)
'random_rotation/rotation_matrix/truedivЧ
'random_rotation/rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2)
'random_rotation/rotation_matrix/sub_5/y┬
%random_rotation/rotation_matrix/sub_5Subrandom_rotation/Cast:y:00random_rotation/rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 2'
%random_rotation/rotation_matrix/sub_5й
%random_rotation/rotation_matrix/Sin_1Sin$random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/Sin_1Ч
'random_rotation/rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2)
'random_rotation/rotation_matrix/sub_6/y─
%random_rotation/rotation_matrix/sub_6Subrandom_rotation/Cast_1:y:00random_rotation/rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 2'
%random_rotation/rotation_matrix/sub_6┘
%random_rotation/rotation_matrix/mul_2Mul)random_rotation/rotation_matrix/Sin_1:y:0)random_rotation/rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/mul_2й
%random_rotation/rotation_matrix/Cos_1Cos$random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/Cos_1Ч
'random_rotation/rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2)
'random_rotation/rotation_matrix/sub_7/y┬
%random_rotation/rotation_matrix/sub_7Subrandom_rotation/Cast:y:00random_rotation/rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 2'
%random_rotation/rotation_matrix/sub_7┘
%random_rotation/rotation_matrix/mul_3Mul)random_rotation/rotation_matrix/Cos_1:y:0)random_rotation/rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/mul_3╫
#random_rotation/rotation_matrix/addAddV2)random_rotation/rotation_matrix/mul_2:z:0)random_rotation/rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:         2%
#random_rotation/rotation_matrix/add╫
%random_rotation/rotation_matrix/sub_8Sub)random_rotation/rotation_matrix/sub_5:z:0'random_rotation/rotation_matrix/add:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/sub_8Я
+random_rotation/rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2-
+random_rotation/rotation_matrix/truediv_1/yЁ
)random_rotation/rotation_matrix/truediv_1RealDiv)random_rotation/rotation_matrix/sub_8:z:04random_rotation/rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:         2+
)random_rotation/rotation_matrix/truediv_1в
%random_rotation/rotation_matrix/ShapeShape$random_rotation/stateful_uniform:z:0*
T0*
_output_shapes
:2'
%random_rotation/rotation_matrix/Shape┤
3random_rotation/rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 25
3random_rotation/rotation_matrix/strided_slice/stack╕
5random_rotation/rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:27
5random_rotation/rotation_matrix/strided_slice/stack_1╕
5random_rotation/rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:27
5random_rotation/rotation_matrix/strided_slice/stack_2в
-random_rotation/rotation_matrix/strided_sliceStridedSlice.random_rotation/rotation_matrix/Shape:output:0<random_rotation/rotation_matrix/strided_slice/stack:output:0>random_rotation/rotation_matrix/strided_slice/stack_1:output:0>random_rotation/rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2/
-random_rotation/rotation_matrix/strided_sliceй
%random_rotation/rotation_matrix/Cos_2Cos$random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/Cos_2┐
5random_rotation/rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        27
5random_rotation/rotation_matrix/strided_slice_1/stack├
7random_rotation/rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation/rotation_matrix/strided_slice_1/stack_1├
7random_rotation/rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      29
7random_rotation/rotation_matrix/strided_slice_1/stack_2╫
/random_rotation/rotation_matrix/strided_slice_1StridedSlice)random_rotation/rotation_matrix/Cos_2:y:0>random_rotation/rotation_matrix/strided_slice_1/stack:output:0@random_rotation/rotation_matrix/strided_slice_1/stack_1:output:0@random_rotation/rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask21
/random_rotation/rotation_matrix/strided_slice_1й
%random_rotation/rotation_matrix/Sin_2Sin$random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/Sin_2┐
5random_rotation/rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        27
5random_rotation/rotation_matrix/strided_slice_2/stack├
7random_rotation/rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation/rotation_matrix/strided_slice_2/stack_1├
7random_rotation/rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      29
7random_rotation/rotation_matrix/strided_slice_2/stack_2╫
/random_rotation/rotation_matrix/strided_slice_2StridedSlice)random_rotation/rotation_matrix/Sin_2:y:0>random_rotation/rotation_matrix/strided_slice_2/stack:output:0@random_rotation/rotation_matrix/strided_slice_2/stack_1:output:0@random_rotation/rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask21
/random_rotation/rotation_matrix/strided_slice_2╜
#random_rotation/rotation_matrix/NegNeg8random_rotation/rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:         2%
#random_rotation/rotation_matrix/Neg┐
5random_rotation/rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        27
5random_rotation/rotation_matrix/strided_slice_3/stack├
7random_rotation/rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation/rotation_matrix/strided_slice_3/stack_1├
7random_rotation/rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      29
7random_rotation/rotation_matrix/strided_slice_3/stack_2┘
/random_rotation/rotation_matrix/strided_slice_3StridedSlice+random_rotation/rotation_matrix/truediv:z:0>random_rotation/rotation_matrix/strided_slice_3/stack:output:0@random_rotation/rotation_matrix/strided_slice_3/stack_1:output:0@random_rotation/rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask21
/random_rotation/rotation_matrix/strided_slice_3й
%random_rotation/rotation_matrix/Sin_3Sin$random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/Sin_3┐
5random_rotation/rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        27
5random_rotation/rotation_matrix/strided_slice_4/stack├
7random_rotation/rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation/rotation_matrix/strided_slice_4/stack_1├
7random_rotation/rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      29
7random_rotation/rotation_matrix/strided_slice_4/stack_2╫
/random_rotation/rotation_matrix/strided_slice_4StridedSlice)random_rotation/rotation_matrix/Sin_3:y:0>random_rotation/rotation_matrix/strided_slice_4/stack:output:0@random_rotation/rotation_matrix/strided_slice_4/stack_1:output:0@random_rotation/rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask21
/random_rotation/rotation_matrix/strided_slice_4й
%random_rotation/rotation_matrix/Cos_3Cos$random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         2'
%random_rotation/rotation_matrix/Cos_3┐
5random_rotation/rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        27
5random_rotation/rotation_matrix/strided_slice_5/stack├
7random_rotation/rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation/rotation_matrix/strided_slice_5/stack_1├
7random_rotation/rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      29
7random_rotation/rotation_matrix/strided_slice_5/stack_2╫
/random_rotation/rotation_matrix/strided_slice_5StridedSlice)random_rotation/rotation_matrix/Cos_3:y:0>random_rotation/rotation_matrix/strided_slice_5/stack:output:0@random_rotation/rotation_matrix/strided_slice_5/stack_1:output:0@random_rotation/rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask21
/random_rotation/rotation_matrix/strided_slice_5┐
5random_rotation/rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        27
5random_rotation/rotation_matrix/strided_slice_6/stack├
7random_rotation/rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        29
7random_rotation/rotation_matrix/strided_slice_6/stack_1├
7random_rotation/rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      29
7random_rotation/rotation_matrix/strided_slice_6/stack_2█
/random_rotation/rotation_matrix/strided_slice_6StridedSlice-random_rotation/rotation_matrix/truediv_1:z:0>random_rotation/rotation_matrix/strided_slice_6/stack:output:0@random_rotation/rotation_matrix/strided_slice_6/stack_1:output:0@random_rotation/rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask21
/random_rotation/rotation_matrix/strided_slice_6Ь
+random_rotation/rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2-
+random_rotation/rotation_matrix/zeros/mul/yь
)random_rotation/rotation_matrix/zeros/mulMul6random_rotation/rotation_matrix/strided_slice:output:04random_rotation/rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2+
)random_rotation/rotation_matrix/zeros/mulЯ
,random_rotation/rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2.
,random_rotation/rotation_matrix/zeros/Less/yч
*random_rotation/rotation_matrix/zeros/LessLess-random_rotation/rotation_matrix/zeros/mul:z:05random_rotation/rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2,
*random_rotation/rotation_matrix/zeros/Lessв
.random_rotation/rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :20
.random_rotation/rotation_matrix/zeros/packed/1Г
,random_rotation/rotation_matrix/zeros/packedPack6random_rotation/rotation_matrix/strided_slice:output:07random_rotation/rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2.
,random_rotation/rotation_matrix/zeros/packedЯ
+random_rotation/rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2-
+random_rotation/rotation_matrix/zeros/Constї
%random_rotation/rotation_matrix/zerosFill5random_rotation/rotation_matrix/zeros/packed:output:04random_rotation/rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:         2'
%random_rotation/rotation_matrix/zerosЬ
+random_rotation/rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2-
+random_rotation/rotation_matrix/concat/axis╚
&random_rotation/rotation_matrix/concatConcatV28random_rotation/rotation_matrix/strided_slice_1:output:0'random_rotation/rotation_matrix/Neg:y:08random_rotation/rotation_matrix/strided_slice_3:output:08random_rotation/rotation_matrix/strided_slice_4:output:08random_rotation/rotation_matrix/strided_slice_5:output:08random_rotation/rotation_matrix/strided_slice_6:output:0.random_rotation/rotation_matrix/zeros:output:04random_rotation/rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2(
&random_rotation/rotation_matrix/concatж
random_rotation/transform/ShapeShape4random_flip/stateless_random_flip_left_right/add:z:0*
T0*
_output_shapes
:2!
random_rotation/transform/Shapeи
-random_rotation/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2/
-random_rotation/transform/strided_slice/stackм
/random_rotation/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:21
/random_rotation/transform/strided_slice/stack_1м
/random_rotation/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:21
/random_rotation/transform/strided_slice/stack_2ъ
'random_rotation/transform/strided_sliceStridedSlice(random_rotation/transform/Shape:output:06random_rotation/transform/strided_slice/stack:output:08random_rotation/transform/strided_slice/stack_1:output:08random_rotation/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2)
'random_rotation/transform/strided_sliceС
$random_rotation/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2&
$random_rotation/transform/fill_value┼
4random_rotation/transform/ImageProjectiveTransformV3ImageProjectiveTransformV34random_flip/stateless_random_flip_left_right/add:z:0/random_rotation/rotation_matrix/concat:output:00random_rotation/transform/strided_slice:output:0-random_rotation/transform/fill_value:output:0*/
_output_shapes
:         &*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR26
4random_rotation/transform/ImageProjectiveTransformV3Я
random_zoom/ShapeShapeIrandom_rotation/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_zoom/ShapeМ
random_zoom/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2!
random_zoom/strided_slice/stackР
!random_zoom/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2#
!random_zoom/strided_slice/stack_1Р
!random_zoom/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2#
!random_zoom/strided_slice/stack_2к
random_zoom/strided_sliceStridedSlicerandom_zoom/Shape:output:0(random_zoom/strided_slice/stack:output:0*random_zoom/strided_slice/stack_1:output:0*random_zoom/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom/strided_sliceЩ
!random_zoom/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:
¤        2#
!random_zoom/strided_slice_1/stackЭ
#random_zoom/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
■        2%
#random_zoom/strided_slice_1/stack_1Ф
#random_zoom/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom/strided_slice_1/stack_2┤
random_zoom/strided_slice_1StridedSlicerandom_zoom/Shape:output:0*random_zoom/strided_slice_1/stack:output:0,random_zoom/strided_slice_1/stack_1:output:0,random_zoom/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom/strided_slice_1В
random_zoom/CastCast$random_zoom/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_zoom/CastЩ
!random_zoom/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2#
!random_zoom/strided_slice_2/stackЭ
#random_zoom/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2%
#random_zoom/strided_slice_2/stack_1Ф
#random_zoom/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2%
#random_zoom/strided_slice_2/stack_2┤
random_zoom/strided_slice_2StridedSlicerandom_zoom/Shape:output:0*random_zoom/strided_slice_2/stack:output:0,random_zoom/strided_slice_2/stack_1:output:0,random_zoom/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
random_zoom/strided_slice_2Ж
random_zoom/Cast_1Cast$random_zoom/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
random_zoom/Cast_1О
$random_zoom/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2&
$random_zoom/stateful_uniform/shape/1╤
"random_zoom/stateful_uniform/shapePack"random_zoom/strided_slice:output:0-random_zoom/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2$
"random_zoom/stateful_uniform/shapeЙ
 random_zoom/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *fff?2"
 random_zoom/stateful_uniform/minЙ
 random_zoom/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *═╠М?2"
 random_zoom/stateful_uniform/maxТ
"random_zoom/stateful_uniform/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2$
"random_zoom/stateful_uniform/Const╔
!random_zoom/stateful_uniform/ProdProd+random_zoom/stateful_uniform/shape:output:0+random_zoom/stateful_uniform/Const:output:0*
T0*
_output_shapes
: 2#
!random_zoom/stateful_uniform/ProdМ
#random_zoom/stateful_uniform/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2%
#random_zoom/stateful_uniform/Cast/xо
#random_zoom/stateful_uniform/Cast_1Cast*random_zoom/stateful_uniform/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2%
#random_zoom/stateful_uniform/Cast_1Х
+random_zoom/stateful_uniform/RngReadAndSkipRngReadAndSkip4random_zoom_stateful_uniform_rngreadandskip_resource,random_zoom/stateful_uniform/Cast/x:output:0'random_zoom/stateful_uniform/Cast_1:y:0*
_output_shapes
:2-
+random_zoom/stateful_uniform/RngReadAndSkipо
0random_zoom/stateful_uniform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 22
0random_zoom/stateful_uniform/strided_slice/stack▓
2random_zoom/stateful_uniform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:24
2random_zoom/stateful_uniform/strided_slice/stack_1▓
2random_zoom/stateful_uniform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:24
2random_zoom/stateful_uniform/strided_slice/stack_2Ц
*random_zoom/stateful_uniform/strided_sliceStridedSlice3random_zoom/stateful_uniform/RngReadAndSkip:value:09random_zoom/stateful_uniform/strided_slice/stack:output:0;random_zoom/stateful_uniform/strided_slice/stack_1:output:0;random_zoom/stateful_uniform/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2,
*random_zoom/stateful_uniform/strided_slice╜
$random_zoom/stateful_uniform/BitcastBitcast3random_zoom/stateful_uniform/strided_slice:output:0*
T0	*
_output_shapes
:*

type02&
$random_zoom/stateful_uniform/Bitcast▓
2random_zoom/stateful_uniform/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:24
2random_zoom/stateful_uniform/strided_slice_1/stack╢
4random_zoom/stateful_uniform/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:26
4random_zoom/stateful_uniform/strided_slice_1/stack_1╢
4random_zoom/stateful_uniform/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:26
4random_zoom/stateful_uniform/strided_slice_1/stack_2О
,random_zoom/stateful_uniform/strided_slice_1StridedSlice3random_zoom/stateful_uniform/RngReadAndSkip:value:0;random_zoom/stateful_uniform/strided_slice_1/stack:output:0=random_zoom/stateful_uniform/strided_slice_1/stack_1:output:0=random_zoom/stateful_uniform/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2.
,random_zoom/stateful_uniform/strided_slice_1├
&random_zoom/stateful_uniform/Bitcast_1Bitcast5random_zoom/stateful_uniform/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02(
&random_zoom/stateful_uniform/Bitcast_1╕
9random_zoom/stateful_uniform/StatelessRandomUniformV2/algConst*
_output_shapes
: *
dtype0*
value	B :2;
9random_zoom/stateful_uniform/StatelessRandomUniformV2/algД
5random_zoom/stateful_uniform/StatelessRandomUniformV2StatelessRandomUniformV2+random_zoom/stateful_uniform/shape:output:0/random_zoom/stateful_uniform/Bitcast_1:output:0-random_zoom/stateful_uniform/Bitcast:output:0Brandom_zoom/stateful_uniform/StatelessRandomUniformV2/alg:output:0*'
_output_shapes
:         27
5random_zoom/stateful_uniform/StatelessRandomUniformV2┬
 random_zoom/stateful_uniform/subSub)random_zoom/stateful_uniform/max:output:0)random_zoom/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2"
 random_zoom/stateful_uniform/subу
 random_zoom/stateful_uniform/mulMul>random_zoom/stateful_uniform/StatelessRandomUniformV2:output:0$random_zoom/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:         2"
 random_zoom/stateful_uniform/mul╚
random_zoom/stateful_uniformAddV2$random_zoom/stateful_uniform/mul:z:0)random_zoom/stateful_uniform/min:output:0*
T0*'
_output_shapes
:         2
random_zoom/stateful_uniformt
random_zoom/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
random_zoom/concat/axis╒
random_zoom/concatConcatV2 random_zoom/stateful_uniform:z:0 random_zoom/stateful_uniform:z:0 random_zoom/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2
random_zoom/concatЙ
random_zoom/zoom_matrix/ShapeShaperandom_zoom/concat:output:0*
T0*
_output_shapes
:2
random_zoom/zoom_matrix/Shapeд
+random_zoom/zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2-
+random_zoom/zoom_matrix/strided_slice/stackи
-random_zoom/zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2/
-random_zoom/zoom_matrix/strided_slice/stack_1и
-random_zoom/zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2/
-random_zoom/zoom_matrix/strided_slice/stack_2Є
%random_zoom/zoom_matrix/strided_sliceStridedSlice&random_zoom/zoom_matrix/Shape:output:04random_zoom/zoom_matrix/strided_slice/stack:output:06random_zoom/zoom_matrix/strided_slice/stack_1:output:06random_zoom/zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2'
%random_zoom/zoom_matrix/strided_sliceГ
random_zoom/zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
random_zoom/zoom_matrix/sub/yв
random_zoom/zoom_matrix/subSubrandom_zoom/Cast_1:y:0&random_zoom/zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 2
random_zoom/zoom_matrix/subЛ
!random_zoom/zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2#
!random_zoom/zoom_matrix/truediv/y╗
random_zoom/zoom_matrix/truedivRealDivrandom_zoom/zoom_matrix/sub:z:0*random_zoom/zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 2!
random_zoom/zoom_matrix/truediv│
-random_zoom/zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2/
-random_zoom/zoom_matrix/strided_slice_1/stack╖
/random_zoom/zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom/zoom_matrix/strided_slice_1/stack_1╖
/random_zoom/zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         21
/random_zoom/zoom_matrix/strided_slice_1/stack_2╣
'random_zoom/zoom_matrix/strided_slice_1StridedSlicerandom_zoom/concat:output:06random_zoom/zoom_matrix/strided_slice_1/stack:output:08random_zoom/zoom_matrix/strided_slice_1/stack_1:output:08random_zoom/zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2)
'random_zoom/zoom_matrix/strided_slice_1З
random_zoom/zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2!
random_zoom/zoom_matrix/sub_1/x╙
random_zoom/zoom_matrix/sub_1Sub(random_zoom/zoom_matrix/sub_1/x:output:00random_zoom/zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:         2
random_zoom/zoom_matrix/sub_1╗
random_zoom/zoom_matrix/mulMul#random_zoom/zoom_matrix/truediv:z:0!random_zoom/zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:         2
random_zoom/zoom_matrix/mulЗ
random_zoom/zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2!
random_zoom/zoom_matrix/sub_2/yж
random_zoom/zoom_matrix/sub_2Subrandom_zoom/Cast:y:0(random_zoom/zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2
random_zoom/zoom_matrix/sub_2П
#random_zoom/zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2%
#random_zoom/zoom_matrix/truediv_1/y├
!random_zoom/zoom_matrix/truediv_1RealDiv!random_zoom/zoom_matrix/sub_2:z:0,random_zoom/zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 2#
!random_zoom/zoom_matrix/truediv_1│
-random_zoom/zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2/
-random_zoom/zoom_matrix/strided_slice_2/stack╖
/random_zoom/zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom/zoom_matrix/strided_slice_2/stack_1╖
/random_zoom/zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         21
/random_zoom/zoom_matrix/strided_slice_2/stack_2╣
'random_zoom/zoom_matrix/strided_slice_2StridedSlicerandom_zoom/concat:output:06random_zoom/zoom_matrix/strided_slice_2/stack:output:08random_zoom/zoom_matrix/strided_slice_2/stack_1:output:08random_zoom/zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2)
'random_zoom/zoom_matrix/strided_slice_2З
random_zoom/zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2!
random_zoom/zoom_matrix/sub_3/x╙
random_zoom/zoom_matrix/sub_3Sub(random_zoom/zoom_matrix/sub_3/x:output:00random_zoom/zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:         2
random_zoom/zoom_matrix/sub_3┴
random_zoom/zoom_matrix/mul_1Mul%random_zoom/zoom_matrix/truediv_1:z:0!random_zoom/zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:         2
random_zoom/zoom_matrix/mul_1│
-random_zoom/zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2/
-random_zoom/zoom_matrix/strided_slice_3/stack╖
/random_zoom/zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom/zoom_matrix/strided_slice_3/stack_1╖
/random_zoom/zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         21
/random_zoom/zoom_matrix/strided_slice_3/stack_2╣
'random_zoom/zoom_matrix/strided_slice_3StridedSlicerandom_zoom/concat:output:06random_zoom/zoom_matrix/strided_slice_3/stack:output:08random_zoom/zoom_matrix/strided_slice_3/stack_1:output:08random_zoom/zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2)
'random_zoom/zoom_matrix/strided_slice_3М
#random_zoom/zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2%
#random_zoom/zoom_matrix/zeros/mul/y╠
!random_zoom/zoom_matrix/zeros/mulMul.random_zoom/zoom_matrix/strided_slice:output:0,random_zoom/zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2#
!random_zoom/zoom_matrix/zeros/mulП
$random_zoom/zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2&
$random_zoom/zoom_matrix/zeros/Less/y╟
"random_zoom/zoom_matrix/zeros/LessLess%random_zoom/zoom_matrix/zeros/mul:z:0-random_zoom/zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2$
"random_zoom/zoom_matrix/zeros/LessТ
&random_zoom/zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2(
&random_zoom/zoom_matrix/zeros/packed/1у
$random_zoom/zoom_matrix/zeros/packedPack.random_zoom/zoom_matrix/strided_slice:output:0/random_zoom/zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2&
$random_zoom/zoom_matrix/zeros/packedП
#random_zoom/zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2%
#random_zoom/zoom_matrix/zeros/Const╒
random_zoom/zoom_matrix/zerosFill-random_zoom/zoom_matrix/zeros/packed:output:0,random_zoom/zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:         2
random_zoom/zoom_matrix/zerosР
%random_zoom/zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2'
%random_zoom/zoom_matrix/zeros_1/mul/y╥
#random_zoom/zoom_matrix/zeros_1/mulMul.random_zoom/zoom_matrix/strided_slice:output:0.random_zoom/zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2%
#random_zoom/zoom_matrix/zeros_1/mulУ
&random_zoom/zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2(
&random_zoom/zoom_matrix/zeros_1/Less/y╧
$random_zoom/zoom_matrix/zeros_1/LessLess'random_zoom/zoom_matrix/zeros_1/mul:z:0/random_zoom/zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2&
$random_zoom/zoom_matrix/zeros_1/LessЦ
(random_zoom/zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2*
(random_zoom/zoom_matrix/zeros_1/packed/1щ
&random_zoom/zoom_matrix/zeros_1/packedPack.random_zoom/zoom_matrix/strided_slice:output:01random_zoom/zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2(
&random_zoom/zoom_matrix/zeros_1/packedУ
%random_zoom/zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%random_zoom/zoom_matrix/zeros_1/Const▌
random_zoom/zoom_matrix/zeros_1Fill/random_zoom/zoom_matrix/zeros_1/packed:output:0.random_zoom/zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:         2!
random_zoom/zoom_matrix/zeros_1│
-random_zoom/zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2/
-random_zoom/zoom_matrix/strided_slice_4/stack╖
/random_zoom/zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           21
/random_zoom/zoom_matrix/strided_slice_4/stack_1╖
/random_zoom/zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         21
/random_zoom/zoom_matrix/strided_slice_4/stack_2╣
'random_zoom/zoom_matrix/strided_slice_4StridedSlicerandom_zoom/concat:output:06random_zoom/zoom_matrix/strided_slice_4/stack:output:08random_zoom/zoom_matrix/strided_slice_4/stack_1:output:08random_zoom/zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2)
'random_zoom/zoom_matrix/strided_slice_4Р
%random_zoom/zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2'
%random_zoom/zoom_matrix/zeros_2/mul/y╥
#random_zoom/zoom_matrix/zeros_2/mulMul.random_zoom/zoom_matrix/strided_slice:output:0.random_zoom/zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2%
#random_zoom/zoom_matrix/zeros_2/mulУ
&random_zoom/zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2(
&random_zoom/zoom_matrix/zeros_2/Less/y╧
$random_zoom/zoom_matrix/zeros_2/LessLess'random_zoom/zoom_matrix/zeros_2/mul:z:0/random_zoom/zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2&
$random_zoom/zoom_matrix/zeros_2/LessЦ
(random_zoom/zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2*
(random_zoom/zoom_matrix/zeros_2/packed/1щ
&random_zoom/zoom_matrix/zeros_2/packedPack.random_zoom/zoom_matrix/strided_slice:output:01random_zoom/zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2(
&random_zoom/zoom_matrix/zeros_2/packedУ
%random_zoom/zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2'
%random_zoom/zoom_matrix/zeros_2/Const▌
random_zoom/zoom_matrix/zeros_2Fill/random_zoom/zoom_matrix/zeros_2/packed:output:0.random_zoom/zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:         2!
random_zoom/zoom_matrix/zeros_2М
#random_zoom/zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2%
#random_zoom/zoom_matrix/concat/axis┘
random_zoom/zoom_matrix/concatConcatV20random_zoom/zoom_matrix/strided_slice_3:output:0&random_zoom/zoom_matrix/zeros:output:0random_zoom/zoom_matrix/mul:z:0(random_zoom/zoom_matrix/zeros_1:output:00random_zoom/zoom_matrix/strided_slice_4:output:0!random_zoom/zoom_matrix/mul_1:z:0(random_zoom/zoom_matrix/zeros_2:output:0,random_zoom/zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2 
random_zoom/zoom_matrix/concat│
random_zoom/transform/ShapeShapeIrandom_rotation/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
random_zoom/transform/Shapeа
)random_zoom/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2+
)random_zoom/transform/strided_slice/stackд
+random_zoom/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2-
+random_zoom/transform/strided_slice/stack_1д
+random_zoom/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2-
+random_zoom/transform/strided_slice/stack_2╥
#random_zoom/transform/strided_sliceStridedSlice$random_zoom/transform/Shape:output:02random_zoom/transform/strided_slice/stack:output:04random_zoom/transform/strided_slice/stack_1:output:04random_zoom/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2%
#random_zoom/transform/strided_sliceЙ
 random_zoom/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 random_zoom/transform/fill_value┬
0random_zoom/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Irandom_rotation/transform/ImageProjectiveTransformV3:transformed_images:0'random_zoom/zoom_matrix/concat:output:0,random_zoom/transform/strided_slice:output:0)random_zoom/transform/fill_value:output:0*/
_output_shapes
:         &*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR22
0random_zoom/transform/ImageProjectiveTransformV3и
IdentityIdentityErandom_zoom/transform/ImageProjectiveTransformV3:transformed_images:0^NoOp*
T0*/
_output_shapes
:         &2

Identityи
NoOpNoOp5^random_flip/stateful_uniform_full_int/RngReadAndSkip\^random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgc^random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter0^random_rotation/stateful_uniform/RngReadAndSkip,^random_zoom/stateful_uniform/RngReadAndSkip*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*4
_input_shapes#
!:         &: : : 2l
4random_flip/stateful_uniform_full_int/RngReadAndSkip4random_flip/stateful_uniform_full_int/RngReadAndSkip2║
[random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg[random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg2╚
brandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterbrandom_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter2b
/random_rotation/stateful_uniform/RngReadAndSkip/random_rotation/stateful_uniform/RngReadAndSkip2Z
+random_zoom/stateful_uniform/RngReadAndSkip+random_zoom/stateful_uniform/RngReadAndSkip:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ю
`
'__inference_dropout_layer_call_fn_78509

inputs
identityИвStatefulPartitionedCall┘
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dropout_layer_call_and_return_conditional_losses_767112
StatefulPartitionedCall|
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:         м2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:         м22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
Ї
Ч
'__inference_dense_1_layer_call_fn_78482

inputs
unknown:
юм
	unknown_0:	м
identityИвStatefulPartitionedCallє
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_1_layer_call_and_return_conditional_losses_764872
StatefulPartitionedCall|
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:         м2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         ю: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:         ю
 
_user_specified_nameinputs
ж
░
__inference_loss_fn_1_78659M
9dense_1_kernel_regularizer_square_readvariableop_resource:
юм
identityИв0dense_1/kernel/Regularizer/Square/ReadVariableOpр
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_1_kernel_regularizer_square_readvariableop_resource* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mull
IdentityIdentity"dense_1/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

IdentityБ
NoOpNoOp1^dense_1/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp
ф
^
B__inference_flatten_layer_call_and_return_conditional_losses_76445

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"    а   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:         а2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:         а2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         P:W S
/
_output_shapes
:         P
 
_user_specified_nameinputs
Ъ
b
F__inference_random_zoom_layer_call_and_return_conditional_losses_78909

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
э
Ф
'__inference_dense_3_layer_call_fn_78573

inputs
unknown:d2
	unknown_0:2
identityИвStatefulPartitionedCallЄ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_765402
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         22

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         d: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         d
 
_user_specified_nameinputs
к
c
D__inference_dropout_1_layer_call_and_return_conditional_losses_76668

inputs
identityИc
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n█╢?2
dropout/Consts
dropout/MulMulinputsdropout/Const:output:0*
T0*'
_output_shapes
:         22
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape┤
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:         2*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЪЩЩ>2
dropout/GreaterEqual/y╛
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:         22
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:         22
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:         22
dropout/Mul_1e
IdentityIdentitydropout/Mul_1:z:0*
T0*'
_output_shapes
:         22

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:         2:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
┤С
┐
F__inference_random_zoom_layer_call_and_return_conditional_losses_75953

inputs6
(stateful_uniform_rngreadandskip_resource:	
identityИвstateful_uniform/RngReadAndSkipD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliceБ
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:
¤        2
strided_slice_1/stackЕ
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ь
strided_slice_1StridedSliceShape:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1^
CastCaststrided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
CastБ
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice_2/stackЕ
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ь
strided_slice_2StridedSliceShape:output:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_2b
Cast_1Caststrided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
Cast_1v
stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2
stateful_uniform/shape/1б
stateful_uniform/shapePackstrided_slice:output:0!stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2
stateful_uniform/shapeq
stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *fff?2
stateful_uniform/minq
stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *═╠М?2
stateful_uniform/maxz
stateful_uniform/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
stateful_uniform/ConstЩ
stateful_uniform/ProdProdstateful_uniform/shape:output:0stateful_uniform/Const:output:0*
T0*
_output_shapes
: 2
stateful_uniform/Prodt
stateful_uniform/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2
stateful_uniform/Cast/xК
stateful_uniform/Cast_1Caststateful_uniform/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
stateful_uniform/Cast_1┘
stateful_uniform/RngReadAndSkipRngReadAndSkip(stateful_uniform_rngreadandskip_resource stateful_uniform/Cast/x:output:0stateful_uniform/Cast_1:y:0*
_output_shapes
:2!
stateful_uniform/RngReadAndSkipЦ
$stateful_uniform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2&
$stateful_uniform/strided_slice/stackЪ
&stateful_uniform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice/stack_1Ъ
&stateful_uniform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice/stack_2╬
stateful_uniform/strided_sliceStridedSlice'stateful_uniform/RngReadAndSkip:value:0-stateful_uniform/strided_slice/stack:output:0/stateful_uniform/strided_slice/stack_1:output:0/stateful_uniform/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2 
stateful_uniform/strided_sliceЩ
stateful_uniform/BitcastBitcast'stateful_uniform/strided_slice:output:0*
T0	*
_output_shapes
:*

type02
stateful_uniform/BitcastЪ
&stateful_uniform/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice_1/stackЮ
(stateful_uniform/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(stateful_uniform/strided_slice_1/stack_1Ю
(stateful_uniform/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(stateful_uniform/strided_slice_1/stack_2╞
 stateful_uniform/strided_slice_1StridedSlice'stateful_uniform/RngReadAndSkip:value:0/stateful_uniform/strided_slice_1/stack:output:01stateful_uniform/strided_slice_1/stack_1:output:01stateful_uniform/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2"
 stateful_uniform/strided_slice_1Я
stateful_uniform/Bitcast_1Bitcast)stateful_uniform/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02
stateful_uniform/Bitcast_1а
-stateful_uniform/StatelessRandomUniformV2/algConst*
_output_shapes
: *
dtype0*
value	B :2/
-stateful_uniform/StatelessRandomUniformV2/alg╝
)stateful_uniform/StatelessRandomUniformV2StatelessRandomUniformV2stateful_uniform/shape:output:0#stateful_uniform/Bitcast_1:output:0!stateful_uniform/Bitcast:output:06stateful_uniform/StatelessRandomUniformV2/alg:output:0*'
_output_shapes
:         2+
)stateful_uniform/StatelessRandomUniformV2Т
stateful_uniform/subSubstateful_uniform/max:output:0stateful_uniform/min:output:0*
T0*
_output_shapes
: 2
stateful_uniform/sub│
stateful_uniform/mulMul2stateful_uniform/StatelessRandomUniformV2:output:0stateful_uniform/sub:z:0*
T0*'
_output_shapes
:         2
stateful_uniform/mulШ
stateful_uniformAddV2stateful_uniform/mul:z:0stateful_uniform/min:output:0*
T0*'
_output_shapes
:         2
stateful_uniform\
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concat/axisЩ
concatConcatV2stateful_uniform:z:0stateful_uniform:z:0concat/axis:output:0*
N*
T0*'
_output_shapes
:         2
concate
zoom_matrix/ShapeShapeconcat:output:0*
T0*
_output_shapes
:2
zoom_matrix/ShapeМ
zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2!
zoom_matrix/strided_slice/stackР
!zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2#
!zoom_matrix/strided_slice/stack_1Р
!zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2#
!zoom_matrix/strided_slice/stack_2к
zoom_matrix/strided_sliceStridedSlicezoom_matrix/Shape:output:0(zoom_matrix/strided_slice/stack:output:0*zoom_matrix/strided_slice/stack_1:output:0*zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
zoom_matrix/strided_slicek
zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
zoom_matrix/sub/yr
zoom_matrix/subSub
Cast_1:y:0zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/subs
zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
zoom_matrix/truediv/yЛ
zoom_matrix/truedivRealDivzoom_matrix/sub:z:0zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/truedivЫ
!zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2#
!zoom_matrix/strided_slice_1/stackЯ
#zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2%
#zoom_matrix/strided_slice_1/stack_1Я
#zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2%
#zoom_matrix/strided_slice_1/stack_2ё
zoom_matrix/strided_slice_1StridedSliceconcat:output:0*zoom_matrix/strided_slice_1/stack:output:0,zoom_matrix/strided_slice_1/stack_1:output:0,zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2
zoom_matrix/strided_slice_1o
zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
zoom_matrix/sub_1/xг
zoom_matrix/sub_1Subzoom_matrix/sub_1/x:output:0$zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/sub_1Л
zoom_matrix/mulMulzoom_matrix/truediv:z:0zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:         2
zoom_matrix/mulo
zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
zoom_matrix/sub_2/yv
zoom_matrix/sub_2SubCast:y:0zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/sub_2w
zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
zoom_matrix/truediv_1/yУ
zoom_matrix/truediv_1RealDivzoom_matrix/sub_2:z:0 zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/truediv_1Ы
!zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2#
!zoom_matrix/strided_slice_2/stackЯ
#zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2%
#zoom_matrix/strided_slice_2/stack_1Я
#zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2%
#zoom_matrix/strided_slice_2/stack_2ё
zoom_matrix/strided_slice_2StridedSliceconcat:output:0*zoom_matrix/strided_slice_2/stack:output:0,zoom_matrix/strided_slice_2/stack_1:output:0,zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2
zoom_matrix/strided_slice_2o
zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
zoom_matrix/sub_3/xг
zoom_matrix/sub_3Subzoom_matrix/sub_3/x:output:0$zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/sub_3С
zoom_matrix/mul_1Mulzoom_matrix/truediv_1:z:0zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:         2
zoom_matrix/mul_1Ы
!zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2#
!zoom_matrix/strided_slice_3/stackЯ
#zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2%
#zoom_matrix/strided_slice_3/stack_1Я
#zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2%
#zoom_matrix/strided_slice_3/stack_2ё
zoom_matrix/strided_slice_3StridedSliceconcat:output:0*zoom_matrix/strided_slice_3/stack:output:0,zoom_matrix/strided_slice_3/stack_1:output:0,zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2
zoom_matrix/strided_slice_3t
zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros/mul/yЬ
zoom_matrix/zeros/mulMul"zoom_matrix/strided_slice:output:0 zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros/mulw
zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zoom_matrix/zeros/Less/yЧ
zoom_matrix/zeros/LessLesszoom_matrix/zeros/mul:z:0!zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros/Lessz
zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros/packed/1│
zoom_matrix/zeros/packedPack"zoom_matrix/strided_slice:output:0#zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zoom_matrix/zeros/packedw
zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zoom_matrix/zeros/Constе
zoom_matrix/zerosFill!zoom_matrix/zeros/packed:output:0 zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/zerosx
zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros_1/mul/yв
zoom_matrix/zeros_1/mulMul"zoom_matrix/strided_slice:output:0"zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros_1/mul{
zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zoom_matrix/zeros_1/Less/yЯ
zoom_matrix/zeros_1/LessLesszoom_matrix/zeros_1/mul:z:0#zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros_1/Less~
zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros_1/packed/1╣
zoom_matrix/zeros_1/packedPack"zoom_matrix/strided_slice:output:0%zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2
zoom_matrix/zeros_1/packed{
zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zoom_matrix/zeros_1/Constн
zoom_matrix/zeros_1Fill#zoom_matrix/zeros_1/packed:output:0"zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/zeros_1Ы
!zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2#
!zoom_matrix/strided_slice_4/stackЯ
#zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2%
#zoom_matrix/strided_slice_4/stack_1Я
#zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2%
#zoom_matrix/strided_slice_4/stack_2ё
zoom_matrix/strided_slice_4StridedSliceconcat:output:0*zoom_matrix/strided_slice_4/stack:output:0,zoom_matrix/strided_slice_4/stack_1:output:0,zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2
zoom_matrix/strided_slice_4x
zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros_2/mul/yв
zoom_matrix/zeros_2/mulMul"zoom_matrix/strided_slice:output:0"zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros_2/mul{
zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zoom_matrix/zeros_2/Less/yЯ
zoom_matrix/zeros_2/LessLesszoom_matrix/zeros_2/mul:z:0#zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros_2/Less~
zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros_2/packed/1╣
zoom_matrix/zeros_2/packedPack"zoom_matrix/strided_slice:output:0%zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2
zoom_matrix/zeros_2/packed{
zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zoom_matrix/zeros_2/Constн
zoom_matrix/zeros_2Fill#zoom_matrix/zeros_2/packed:output:0"zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/zeros_2t
zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/concat/axisс
zoom_matrix/concatConcatV2$zoom_matrix/strided_slice_3:output:0zoom_matrix/zeros:output:0zoom_matrix/mul:z:0zoom_matrix/zeros_1:output:0$zoom_matrix/strided_slice_4:output:0zoom_matrix/mul_1:z:0zoom_matrix/zeros_2:output:0 zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2
zoom_matrix/concatX
transform/ShapeShapeinputs*
T0*
_output_shapes
:2
transform/ShapeИ
transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2
transform/strided_slice/stackМ
transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2!
transform/strided_slice/stack_1М
transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2!
transform/strided_slice/stack_2К
transform/strided_sliceStridedSlicetransform/Shape:output:0&transform/strided_slice/stack:output:0(transform/strided_slice/stack_1:output:0(transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2
transform/strided_sliceq
transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2
transform/fill_value├
$transform/ImageProjectiveTransformV3ImageProjectiveTransformV3inputszoom_matrix/concat:output:0 transform/strided_slice:output:0transform/fill_value:output:0*/
_output_shapes
:         &*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2&
$transform/ImageProjectiveTransformV3Ь
IdentityIdentity9transform/ImageProjectiveTransformV3:transformed_images:0^NoOp*
T0*/
_output_shapes
:         &2

Identityp
NoOpNoOp ^stateful_uniform/RngReadAndSkip*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 2B
stateful_uniform/RngReadAndSkipstateful_uniform/RngReadAndSkip:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Щ
a
E__inference_sequential_layer_call_and_return_conditional_losses_77945

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
щ
№
C__inference_conv2d_4_layer_call_and_return_conditional_losses_78404

inputs8
conv2d_readvariableop_resource:@P-
biasadd_readvariableop_resource:P
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:@P*
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P*
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:P*
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         P2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         P2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         @: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:         @
 
_user_specified_nameinputs
к
f
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_76255

inputs
identityн
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4                                    *
ksize
*
paddingVALID*
strides
2	
MaxPoolЗ
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
╨f
№
F__inference_random_flip_layer_call_and_return_conditional_losses_76159

inputs?
1stateful_uniform_full_int_rngreadandskip_resource:	
identityИв(stateful_uniform_full_int/RngReadAndSkipвOstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgвVstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterМ
stateful_uniform_full_int/shapeConst*
_output_shapes
:*
dtype0*
valueB:2!
stateful_uniform_full_int/shapeМ
stateful_uniform_full_int/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
stateful_uniform_full_int/Const╜
stateful_uniform_full_int/ProdProd(stateful_uniform_full_int/shape:output:0(stateful_uniform_full_int/Const:output:0*
T0*
_output_shapes
: 2 
stateful_uniform_full_int/ProdЖ
 stateful_uniform_full_int/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2"
 stateful_uniform_full_int/Cast/xе
 stateful_uniform_full_int/Cast_1Cast'stateful_uniform_full_int/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2"
 stateful_uniform_full_int/Cast_1Ж
(stateful_uniform_full_int/RngReadAndSkipRngReadAndSkip1stateful_uniform_full_int_rngreadandskip_resource)stateful_uniform_full_int/Cast/x:output:0$stateful_uniform_full_int/Cast_1:y:0*
_output_shapes
:2*
(stateful_uniform_full_int/RngReadAndSkipи
-stateful_uniform_full_int/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2/
-stateful_uniform_full_int/strided_slice/stackм
/stateful_uniform_full_int/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:21
/stateful_uniform_full_int/strided_slice/stack_1м
/stateful_uniform_full_int/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:21
/stateful_uniform_full_int/strided_slice/stack_2Д
'stateful_uniform_full_int/strided_sliceStridedSlice0stateful_uniform_full_int/RngReadAndSkip:value:06stateful_uniform_full_int/strided_slice/stack:output:08stateful_uniform_full_int/strided_slice/stack_1:output:08stateful_uniform_full_int/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2)
'stateful_uniform_full_int/strided_slice┤
!stateful_uniform_full_int/BitcastBitcast0stateful_uniform_full_int/strided_slice:output:0*
T0	*
_output_shapes
:*

type02#
!stateful_uniform_full_int/Bitcastм
/stateful_uniform_full_int/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:21
/stateful_uniform_full_int/strided_slice_1/stack░
1stateful_uniform_full_int/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:23
1stateful_uniform_full_int/strided_slice_1/stack_1░
1stateful_uniform_full_int/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:23
1stateful_uniform_full_int/strided_slice_1/stack_2№
)stateful_uniform_full_int/strided_slice_1StridedSlice0stateful_uniform_full_int/RngReadAndSkip:value:08stateful_uniform_full_int/strided_slice_1/stack:output:0:stateful_uniform_full_int/strided_slice_1/stack_1:output:0:stateful_uniform_full_int/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2+
)stateful_uniform_full_int/strided_slice_1║
#stateful_uniform_full_int/Bitcast_1Bitcast2stateful_uniform_full_int/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02%
#stateful_uniform_full_int/Bitcast_1А
stateful_uniform_full_int/algConst*
_output_shapes
: *
dtype0*
value	B :2
stateful_uniform_full_int/algо
stateful_uniform_full_intStatelessRandomUniformFullIntV2(stateful_uniform_full_int/shape:output:0,stateful_uniform_full_int/Bitcast_1:output:0*stateful_uniform_full_int/Bitcast:output:0&stateful_uniform_full_int/alg:output:0*
_output_shapes
:*
dtype0	2
stateful_uniform_full_intb

zeros_likeConst*
_output_shapes
:*
dtype0	*
valueB	R 2

zeros_likeБ
stackPack"stateful_uniform_full_int:output:0zeros_like:output:0*
N*
T0	*
_output_shapes

:2
stack{
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        2
strided_slice/stack
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       2
strided_slice/stack_1
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2
strided_slice/stack_2И
strided_sliceStridedSlicestack:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask*
end_mask*
shrink_axis_mask2
strided_slice╙
3stateless_random_flip_left_right/control_dependencyIdentityinputs*
T0*
_class
loc:@inputs*/
_output_shapes
:         &25
3stateless_random_flip_left_right/control_dependency╝
&stateless_random_flip_left_right/ShapeShape<stateless_random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:2(
&stateless_random_flip_left_right/Shape╢
4stateless_random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 26
4stateless_random_flip_left_right/strided_slice/stack║
6stateless_random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:28
6stateless_random_flip_left_right/strided_slice/stack_1║
6stateless_random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:28
6stateless_random_flip_left_right/strided_slice/stack_2и
.stateless_random_flip_left_right/strided_sliceStridedSlice/stateless_random_flip_left_right/Shape:output:0=stateless_random_flip_left_right/strided_slice/stack:output:0?stateless_random_flip_left_right/strided_slice/stack_1:output:0?stateless_random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask20
.stateless_random_flip_left_right/strided_sliceё
?stateless_random_flip_left_right/stateless_random_uniform/shapePack7stateless_random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2A
?stateless_random_flip_left_right/stateless_random_uniform/shape├
=stateless_random_flip_left_right/stateless_random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2?
=stateless_random_flip_left_right/stateless_random_uniform/min├
=stateless_random_flip_left_right/stateless_random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2?
=stateless_random_flip_left_right/stateless_random_uniform/maxК
Vstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterStatelessRandomGetKeyCounterstrided_slice:output:0* 
_output_shapes
::2X
Vstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterм
Ostateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgStatelessRandomGetAlgW^stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter*
_output_shapes
: 2Q
Ostateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg╩
Rstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2StatelessRandomUniformV2Hstateless_random_flip_left_right/stateless_random_uniform/shape:output:0\stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:key:0`stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:counter:0Ustateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg:alg:0*#
_output_shapes
:         2T
Rstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2╢
=stateless_random_flip_left_right/stateless_random_uniform/subSubFstateless_random_flip_left_right/stateless_random_uniform/max:output:0Fstateless_random_flip_left_right/stateless_random_uniform/min:output:0*
T0*
_output_shapes
: 2?
=stateless_random_flip_left_right/stateless_random_uniform/sub╙
=stateless_random_flip_left_right/stateless_random_uniform/mulMul[stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2:output:0Astateless_random_flip_left_right/stateless_random_uniform/sub:z:0*
T0*#
_output_shapes
:         2?
=stateless_random_flip_left_right/stateless_random_uniform/mul╕
9stateless_random_flip_left_right/stateless_random_uniformAddV2Astateless_random_flip_left_right/stateless_random_uniform/mul:z:0Fstateless_random_flip_left_right/stateless_random_uniform/min:output:0*
T0*#
_output_shapes
:         2;
9stateless_random_flip_left_right/stateless_random_uniformж
0stateless_random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :22
0stateless_random_flip_left_right/Reshape/shape/1ж
0stateless_random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :22
0stateless_random_flip_left_right/Reshape/shape/2ж
0stateless_random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :22
0stateless_random_flip_left_right/Reshape/shape/3А
.stateless_random_flip_left_right/Reshape/shapePack7stateless_random_flip_left_right/strided_slice:output:09stateless_random_flip_left_right/Reshape/shape/1:output:09stateless_random_flip_left_right/Reshape/shape/2:output:09stateless_random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:20
.stateless_random_flip_left_right/Reshape/shapeС
(stateless_random_flip_left_right/ReshapeReshape=stateless_random_flip_left_right/stateless_random_uniform:z:07stateless_random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:         2*
(stateless_random_flip_left_right/Reshape╞
&stateless_random_flip_left_right/RoundRound1stateless_random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:         2(
&stateless_random_flip_left_right/Roundм
/stateless_random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:21
/stateless_random_flip_left_right/ReverseV2/axisЧ
*stateless_random_flip_left_right/ReverseV2	ReverseV2<stateless_random_flip_left_right/control_dependency:output:08stateless_random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:         &2,
*stateless_random_flip_left_right/ReverseV2ю
$stateless_random_flip_left_right/mulMul*stateless_random_flip_left_right/Round:y:03stateless_random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:         &2&
$stateless_random_flip_left_right/mulХ
&stateless_random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2(
&stateless_random_flip_left_right/sub/xъ
$stateless_random_flip_left_right/subSub/stateless_random_flip_left_right/sub/x:output:0*stateless_random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:         2&
$stateless_random_flip_left_right/sub∙
&stateless_random_flip_left_right/mul_1Mul(stateless_random_flip_left_right/sub:z:0<stateless_random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:         &2(
&stateless_random_flip_left_right/mul_1х
$stateless_random_flip_left_right/addAddV2(stateless_random_flip_left_right/mul:z:0*stateless_random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:         &2&
$stateless_random_flip_left_right/addЛ
IdentityIdentity(stateless_random_flip_left_right/add:z:0^NoOp*
T0*/
_output_shapes
:         &2

Identityд
NoOpNoOp)^stateful_uniform_full_int/RngReadAndSkipP^stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgW^stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 2T
(stateful_uniform_full_int/RngReadAndSkip(stateful_uniform_full_int/RngReadAndSkip2в
Ostateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgOstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg2░
Vstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterVstateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
иа
├
J__inference_random_rotation_layer_call_and_return_conditional_losses_78893

inputs6
(stateful_uniform_rngreadandskip_resource:	
identityИвstateful_uniform/RngReadAndSkipD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliceБ
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:
¤        2
strided_slice_1/stackЕ
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ь
strided_slice_1StridedSliceShape:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1^
CastCaststrided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
CastБ
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice_2/stackЕ
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ь
strided_slice_2StridedSliceShape:output:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_2b
Cast_1Caststrided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
Cast_1~
stateful_uniform/shapePackstrided_slice:output:0*
N*
T0*
_output_shapes
:2
stateful_uniform/shapeq
stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *█╔└2
stateful_uniform/minq
stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *█╔@2
stateful_uniform/maxz
stateful_uniform/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
stateful_uniform/ConstЩ
stateful_uniform/ProdProdstateful_uniform/shape:output:0stateful_uniform/Const:output:0*
T0*
_output_shapes
: 2
stateful_uniform/Prodt
stateful_uniform/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2
stateful_uniform/Cast/xК
stateful_uniform/Cast_1Caststateful_uniform/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
stateful_uniform/Cast_1┘
stateful_uniform/RngReadAndSkipRngReadAndSkip(stateful_uniform_rngreadandskip_resource stateful_uniform/Cast/x:output:0stateful_uniform/Cast_1:y:0*
_output_shapes
:2!
stateful_uniform/RngReadAndSkipЦ
$stateful_uniform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2&
$stateful_uniform/strided_slice/stackЪ
&stateful_uniform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice/stack_1Ъ
&stateful_uniform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice/stack_2╬
stateful_uniform/strided_sliceStridedSlice'stateful_uniform/RngReadAndSkip:value:0-stateful_uniform/strided_slice/stack:output:0/stateful_uniform/strided_slice/stack_1:output:0/stateful_uniform/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2 
stateful_uniform/strided_sliceЩ
stateful_uniform/BitcastBitcast'stateful_uniform/strided_slice:output:0*
T0	*
_output_shapes
:*

type02
stateful_uniform/BitcastЪ
&stateful_uniform/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice_1/stackЮ
(stateful_uniform/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(stateful_uniform/strided_slice_1/stack_1Ю
(stateful_uniform/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(stateful_uniform/strided_slice_1/stack_2╞
 stateful_uniform/strided_slice_1StridedSlice'stateful_uniform/RngReadAndSkip:value:0/stateful_uniform/strided_slice_1/stack:output:01stateful_uniform/strided_slice_1/stack_1:output:01stateful_uniform/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2"
 stateful_uniform/strided_slice_1Я
stateful_uniform/Bitcast_1Bitcast)stateful_uniform/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02
stateful_uniform/Bitcast_1а
-stateful_uniform/StatelessRandomUniformV2/algConst*
_output_shapes
: *
dtype0*
value	B :2/
-stateful_uniform/StatelessRandomUniformV2/alg╕
)stateful_uniform/StatelessRandomUniformV2StatelessRandomUniformV2stateful_uniform/shape:output:0#stateful_uniform/Bitcast_1:output:0!stateful_uniform/Bitcast:output:06stateful_uniform/StatelessRandomUniformV2/alg:output:0*#
_output_shapes
:         2+
)stateful_uniform/StatelessRandomUniformV2Т
stateful_uniform/subSubstateful_uniform/max:output:0stateful_uniform/min:output:0*
T0*
_output_shapes
: 2
stateful_uniform/subп
stateful_uniform/mulMul2stateful_uniform/StatelessRandomUniformV2:output:0stateful_uniform/sub:z:0*
T0*#
_output_shapes
:         2
stateful_uniform/mulФ
stateful_uniformAddV2stateful_uniform/mul:z:0stateful_uniform/min:output:0*
T0*#
_output_shapes
:         2
stateful_uniforms
rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub/y~
rotation_matrix/subSub
Cast_1:y:0rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/subu
rotation_matrix/CosCosstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Cosw
rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_1/yД
rotation_matrix/sub_1Sub
Cast_1:y:0 rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_1У
rotation_matrix/mulMulrotation_matrix/Cos:y:0rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/mulu
rotation_matrix/SinSinstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Sinw
rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_2/yВ
rotation_matrix/sub_2SubCast:y:0 rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_2Ч
rotation_matrix/mul_1Mulrotation_matrix/Sin:y:0rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/mul_1Ч
rotation_matrix/sub_3Subrotation_matrix/mul:z:0rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/sub_3Ч
rotation_matrix/sub_4Subrotation_matrix/sub:z:0rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/sub_4{
rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
rotation_matrix/truediv/yк
rotation_matrix/truedivRealDivrotation_matrix/sub_4:z:0"rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:         2
rotation_matrix/truedivw
rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_5/yВ
rotation_matrix/sub_5SubCast:y:0 rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_5y
rotation_matrix/Sin_1Sinstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Sin_1w
rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_6/yД
rotation_matrix/sub_6Sub
Cast_1:y:0 rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_6Щ
rotation_matrix/mul_2Mulrotation_matrix/Sin_1:y:0rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/mul_2y
rotation_matrix/Cos_1Cosstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Cos_1w
rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
rotation_matrix/sub_7/yВ
rotation_matrix/sub_7SubCast:y:0 rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/sub_7Щ
rotation_matrix/mul_3Mulrotation_matrix/Cos_1:y:0rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/mul_3Ч
rotation_matrix/addAddV2rotation_matrix/mul_2:z:0rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/addЧ
rotation_matrix/sub_8Subrotation_matrix/sub_5:z:0rotation_matrix/add:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/sub_8
rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
rotation_matrix/truediv_1/y░
rotation_matrix/truediv_1RealDivrotation_matrix/sub_8:z:0$rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:         2
rotation_matrix/truediv_1r
rotation_matrix/ShapeShapestateful_uniform:z:0*
T0*
_output_shapes
:2
rotation_matrix/ShapeФ
#rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2%
#rotation_matrix/strided_slice/stackШ
%rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2'
%rotation_matrix/strided_slice/stack_1Ш
%rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%rotation_matrix/strided_slice/stack_2┬
rotation_matrix/strided_sliceStridedSlicerotation_matrix/Shape:output:0,rotation_matrix/strided_slice/stack:output:0.rotation_matrix/strided_slice/stack_1:output:0.rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
rotation_matrix/strided_slicey
rotation_matrix/Cos_2Cosstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Cos_2Я
%rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_1/stackг
'rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_1/stack_1г
'rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_1/stack_2ў
rotation_matrix/strided_slice_1StridedSlicerotation_matrix/Cos_2:y:0.rotation_matrix/strided_slice_1/stack:output:00rotation_matrix/strided_slice_1/stack_1:output:00rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_1y
rotation_matrix/Sin_2Sinstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Sin_2Я
%rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_2/stackг
'rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_2/stack_1г
'rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_2/stack_2ў
rotation_matrix/strided_slice_2StridedSlicerotation_matrix/Sin_2:y:0.rotation_matrix/strided_slice_2/stack:output:00rotation_matrix/strided_slice_2/stack_1:output:00rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_2Н
rotation_matrix/NegNeg(rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:         2
rotation_matrix/NegЯ
%rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_3/stackг
'rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_3/stack_1г
'rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_3/stack_2∙
rotation_matrix/strided_slice_3StridedSlicerotation_matrix/truediv:z:0.rotation_matrix/strided_slice_3/stack:output:00rotation_matrix/strided_slice_3/stack_1:output:00rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_3y
rotation_matrix/Sin_3Sinstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Sin_3Я
%rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_4/stackг
'rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_4/stack_1г
'rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_4/stack_2ў
rotation_matrix/strided_slice_4StridedSlicerotation_matrix/Sin_3:y:0.rotation_matrix/strided_slice_4/stack:output:00rotation_matrix/strided_slice_4/stack_1:output:00rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_4y
rotation_matrix/Cos_3Cosstateful_uniform:z:0*
T0*#
_output_shapes
:         2
rotation_matrix/Cos_3Я
%rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_5/stackг
'rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_5/stack_1г
'rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_5/stack_2ў
rotation_matrix/strided_slice_5StridedSlicerotation_matrix/Cos_3:y:0.rotation_matrix/strided_slice_5/stack:output:00rotation_matrix/strided_slice_5/stack_1:output:00rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_5Я
%rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        2'
%rotation_matrix/strided_slice_6/stackг
'rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2)
'rotation_matrix/strided_slice_6/stack_1г
'rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2)
'rotation_matrix/strided_slice_6/stack_2√
rotation_matrix/strided_slice_6StridedSlicerotation_matrix/truediv_1:z:0.rotation_matrix/strided_slice_6/stack:output:00rotation_matrix/strided_slice_6/stack_1:output:00rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2!
rotation_matrix/strided_slice_6|
rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
rotation_matrix/zeros/mul/yм
rotation_matrix/zeros/mulMul&rotation_matrix/strided_slice:output:0$rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/zeros/mul
rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
rotation_matrix/zeros/Less/yз
rotation_matrix/zeros/LessLessrotation_matrix/zeros/mul:z:0%rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
rotation_matrix/zeros/LessВ
rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2 
rotation_matrix/zeros/packed/1├
rotation_matrix/zeros/packedPack&rotation_matrix/strided_slice:output:0'rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
rotation_matrix/zeros/packed
rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rotation_matrix/zeros/Const╡
rotation_matrix/zerosFill%rotation_matrix/zeros/packed:output:0$rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:         2
rotation_matrix/zeros|
rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
rotation_matrix/concat/axisи
rotation_matrix/concatConcatV2(rotation_matrix/strided_slice_1:output:0rotation_matrix/Neg:y:0(rotation_matrix/strided_slice_3:output:0(rotation_matrix/strided_slice_4:output:0(rotation_matrix/strided_slice_5:output:0(rotation_matrix/strided_slice_6:output:0rotation_matrix/zeros:output:0$rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2
rotation_matrix/concatX
transform/ShapeShapeinputs*
T0*
_output_shapes
:2
transform/ShapeИ
transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2
transform/strided_slice/stackМ
transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2!
transform/strided_slice/stack_1М
transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2!
transform/strided_slice/stack_2К
transform/strided_sliceStridedSlicetransform/Shape:output:0&transform/strided_slice/stack:output:0(transform/strided_slice/stack_1:output:0(transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2
transform/strided_sliceq
transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2
transform/fill_value╟
$transform/ImageProjectiveTransformV3ImageProjectiveTransformV3inputsrotation_matrix/concat:output:0 transform/strided_slice:output:0transform/fill_value:output:0*/
_output_shapes
:         &*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2&
$transform/ImageProjectiveTransformV3Ь
IdentityIdentity9transform/ImageProjectiveTransformV3:transformed_images:0^NoOp*
T0*/
_output_shapes
:         &2

Identityp
NoOpNoOp ^stateful_uniform/RngReadAndSkip*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 2B
stateful_uniform/RngReadAndSkipstateful_uniform/RngReadAndSkip:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
╒
K
/__inference_max_pooling2d_3_layer_call_fn_78409

inputs
identityы
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4                                    * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_762992
PartitionedCallП
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
╒
П
,__inference_sequential_1_layer_call_fn_77056
sequential_input
unknown:	
	unknown_0:	
	unknown_1:	#
	unknown_2:
	unknown_3:#
	unknown_4:
	unknown_5:#
	unknown_6: 
	unknown_7: #
	unknown_8: @
	unknown_9:@$

unknown_10:@P

unknown_11:P

unknown_12:
аю

unknown_13:	ю

unknown_14:
юм

unknown_15:	м

unknown_16:	мd

unknown_17:d

unknown_18:d2

unknown_19:2

unknown_20:2

unknown_21:
identityИвStatefulPartitionedCallЫ
StatefulPartitionedCallStatefulPartitionedCallsequential_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21*#
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8В *P
fKRI
G__inference_sequential_1_layer_call_and_return_conditional_losses_769562
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*\
_input_shapesK
I:         &: : : : : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:a ]
/
_output_shapes
:         &
*
_user_specified_namesequential_input
щ
№
C__inference_conv2d_4_layer_call_and_return_conditional_losses_76427

inputs8
conv2d_readvariableop_resource:@P-
biasadd_readvariableop_resource:P
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:@P*
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P*
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:P*
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         P2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         P2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         @: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:         @
 
_user_specified_nameinputs
┌Л
╞
__inference__traced_save_79277
file_prefix,
(savev2_conv2d_kernel_read_readvariableop*
&savev2_conv2d_bias_read_readvariableop.
*savev2_conv2d_1_kernel_read_readvariableop,
(savev2_conv2d_1_bias_read_readvariableop.
*savev2_conv2d_2_kernel_read_readvariableop,
(savev2_conv2d_2_bias_read_readvariableop.
*savev2_conv2d_3_kernel_read_readvariableop,
(savev2_conv2d_3_bias_read_readvariableop.
*savev2_conv2d_4_kernel_read_readvariableop,
(savev2_conv2d_4_bias_read_readvariableop+
'savev2_dense_kernel_read_readvariableop)
%savev2_dense_bias_read_readvariableop-
)savev2_dense_1_kernel_read_readvariableop+
'savev2_dense_1_bias_read_readvariableop-
)savev2_dense_2_kernel_read_readvariableop+
'savev2_dense_2_bias_read_readvariableop-
)savev2_dense_3_kernel_read_readvariableop+
'savev2_dense_3_bias_read_readvariableop,
(savev2_output_kernel_read_readvariableop*
&savev2_output_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop'
#savev2_variable_read_readvariableop	)
%savev2_variable_1_read_readvariableop	)
%savev2_variable_2_read_readvariableop	$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop3
/savev2_adam_conv2d_kernel_m_read_readvariableop1
-savev2_adam_conv2d_bias_m_read_readvariableop5
1savev2_adam_conv2d_1_kernel_m_read_readvariableop3
/savev2_adam_conv2d_1_bias_m_read_readvariableop5
1savev2_adam_conv2d_2_kernel_m_read_readvariableop3
/savev2_adam_conv2d_2_bias_m_read_readvariableop5
1savev2_adam_conv2d_3_kernel_m_read_readvariableop3
/savev2_adam_conv2d_3_bias_m_read_readvariableop5
1savev2_adam_conv2d_4_kernel_m_read_readvariableop3
/savev2_adam_conv2d_4_bias_m_read_readvariableop2
.savev2_adam_dense_kernel_m_read_readvariableop0
,savev2_adam_dense_bias_m_read_readvariableop4
0savev2_adam_dense_1_kernel_m_read_readvariableop2
.savev2_adam_dense_1_bias_m_read_readvariableop4
0savev2_adam_dense_2_kernel_m_read_readvariableop2
.savev2_adam_dense_2_bias_m_read_readvariableop4
0savev2_adam_dense_3_kernel_m_read_readvariableop2
.savev2_adam_dense_3_bias_m_read_readvariableop3
/savev2_adam_output_kernel_m_read_readvariableop1
-savev2_adam_output_bias_m_read_readvariableop3
/savev2_adam_conv2d_kernel_v_read_readvariableop1
-savev2_adam_conv2d_bias_v_read_readvariableop5
1savev2_adam_conv2d_1_kernel_v_read_readvariableop3
/savev2_adam_conv2d_1_bias_v_read_readvariableop5
1savev2_adam_conv2d_2_kernel_v_read_readvariableop3
/savev2_adam_conv2d_2_bias_v_read_readvariableop5
1savev2_adam_conv2d_3_kernel_v_read_readvariableop3
/savev2_adam_conv2d_3_bias_v_read_readvariableop5
1savev2_adam_conv2d_4_kernel_v_read_readvariableop3
/savev2_adam_conv2d_4_bias_v_read_readvariableop2
.savev2_adam_dense_kernel_v_read_readvariableop0
,savev2_adam_dense_bias_v_read_readvariableop4
0savev2_adam_dense_1_kernel_v_read_readvariableop2
.savev2_adam_dense_1_bias_v_read_readvariableop4
0savev2_adam_dense_2_kernel_v_read_readvariableop2
.savev2_adam_dense_2_bias_v_read_readvariableop4
0savev2_adam_dense_3_kernel_v_read_readvariableop2
.savev2_adam_dense_3_bias_v_read_readvariableop3
/savev2_adam_output_kernel_v_read_readvariableop1
-savev2_adam_output_bias_v_read_readvariableop
savev2_const

identity_1ИвMergeV2CheckpointsП
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
Const_1Л
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
ShardedFilename/shardж
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename╓(
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:I*
dtype0*ш'
value▐'B█'IB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-9/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-9/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB:layer-0/layer-0/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEB:layer-0/layer-1/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEB:layer-0/layer-2/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-9/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-9/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-9/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-9/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_namesЭ
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:I*
dtype0*з
valueЭBЪIB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slices╜
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0(savev2_conv2d_kernel_read_readvariableop&savev2_conv2d_bias_read_readvariableop*savev2_conv2d_1_kernel_read_readvariableop(savev2_conv2d_1_bias_read_readvariableop*savev2_conv2d_2_kernel_read_readvariableop(savev2_conv2d_2_bias_read_readvariableop*savev2_conv2d_3_kernel_read_readvariableop(savev2_conv2d_3_bias_read_readvariableop*savev2_conv2d_4_kernel_read_readvariableop(savev2_conv2d_4_bias_read_readvariableop'savev2_dense_kernel_read_readvariableop%savev2_dense_bias_read_readvariableop)savev2_dense_1_kernel_read_readvariableop'savev2_dense_1_bias_read_readvariableop)savev2_dense_2_kernel_read_readvariableop'savev2_dense_2_bias_read_readvariableop)savev2_dense_3_kernel_read_readvariableop'savev2_dense_3_bias_read_readvariableop(savev2_output_kernel_read_readvariableop&savev2_output_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop#savev2_variable_read_readvariableop%savev2_variable_1_read_readvariableop%savev2_variable_2_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop/savev2_adam_conv2d_kernel_m_read_readvariableop-savev2_adam_conv2d_bias_m_read_readvariableop1savev2_adam_conv2d_1_kernel_m_read_readvariableop/savev2_adam_conv2d_1_bias_m_read_readvariableop1savev2_adam_conv2d_2_kernel_m_read_readvariableop/savev2_adam_conv2d_2_bias_m_read_readvariableop1savev2_adam_conv2d_3_kernel_m_read_readvariableop/savev2_adam_conv2d_3_bias_m_read_readvariableop1savev2_adam_conv2d_4_kernel_m_read_readvariableop/savev2_adam_conv2d_4_bias_m_read_readvariableop.savev2_adam_dense_kernel_m_read_readvariableop,savev2_adam_dense_bias_m_read_readvariableop0savev2_adam_dense_1_kernel_m_read_readvariableop.savev2_adam_dense_1_bias_m_read_readvariableop0savev2_adam_dense_2_kernel_m_read_readvariableop.savev2_adam_dense_2_bias_m_read_readvariableop0savev2_adam_dense_3_kernel_m_read_readvariableop.savev2_adam_dense_3_bias_m_read_readvariableop/savev2_adam_output_kernel_m_read_readvariableop-savev2_adam_output_bias_m_read_readvariableop/savev2_adam_conv2d_kernel_v_read_readvariableop-savev2_adam_conv2d_bias_v_read_readvariableop1savev2_adam_conv2d_1_kernel_v_read_readvariableop/savev2_adam_conv2d_1_bias_v_read_readvariableop1savev2_adam_conv2d_2_kernel_v_read_readvariableop/savev2_adam_conv2d_2_bias_v_read_readvariableop1savev2_adam_conv2d_3_kernel_v_read_readvariableop/savev2_adam_conv2d_3_bias_v_read_readvariableop1savev2_adam_conv2d_4_kernel_v_read_readvariableop/savev2_adam_conv2d_4_bias_v_read_readvariableop.savev2_adam_dense_kernel_v_read_readvariableop,savev2_adam_dense_bias_v_read_readvariableop0savev2_adam_dense_1_kernel_v_read_readvariableop.savev2_adam_dense_1_bias_v_read_readvariableop0savev2_adam_dense_2_kernel_v_read_readvariableop.savev2_adam_dense_2_bias_v_read_readvariableop0savev2_adam_dense_3_kernel_v_read_readvariableop.savev2_adam_dense_3_bias_v_read_readvariableop/savev2_adam_output_kernel_v_read_readvariableop-savev2_adam_output_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *W
dtypesM
K2I				2
SaveV2║
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixesб
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identity_

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: 2

Identity_1c
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"!

identity_1Identity_1:output:0*к
_input_shapesШ
Х: ::::: : : @:@:@P:P:
аю:ю:
юм:м:	мd:d:d2:2:2:: : : : : :::: : : : ::::: : : @:@:@P:P:
аю:ю:
юм:м:	мd:d:d2:2:2:::::: : : @:@:@P:P:
аю:ю:
юм:м:	мd:d:d2:2:2:: 2(
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
: :,(
&
_output_shapes
: @: 

_output_shapes
:@:,	(
&
_output_shapes
:@P: 


_output_shapes
:P:&"
 
_output_shapes
:
аю:!

_output_shapes	
:ю:&"
 
_output_shapes
:
юм:!

_output_shapes	
:м:%!

_output_shapes
:	мd: 

_output_shapes
:d:$ 

_output_shapes

:d2: 

_output_shapes
:2:$ 

_output_shapes

:2: 

_output_shapes
::

_output_shapes
: :
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
: : 

_output_shapes
:: 

_output_shapes
:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: : 

_output_shapes
: :,!(
&
_output_shapes
:: "

_output_shapes
::,#(
&
_output_shapes
:: $

_output_shapes
::,%(
&
_output_shapes
: : &

_output_shapes
: :,'(
&
_output_shapes
: @: (

_output_shapes
:@:,)(
&
_output_shapes
:@P: *

_output_shapes
:P:&+"
 
_output_shapes
:
аю:!,

_output_shapes	
:ю:&-"
 
_output_shapes
:
юм:!.

_output_shapes	
:м:%/!

_output_shapes
:	мd: 0

_output_shapes
:d:$1 

_output_shapes

:d2: 2

_output_shapes
:2:$3 

_output_shapes

:2: 4

_output_shapes
::,5(
&
_output_shapes
:: 6

_output_shapes
::,7(
&
_output_shapes
:: 8

_output_shapes
::,9(
&
_output_shapes
: : :

_output_shapes
: :,;(
&
_output_shapes
: @: <

_output_shapes
:@:,=(
&
_output_shapes
:@P: >

_output_shapes
:P:&?"
 
_output_shapes
:
аю:!@

_output_shapes	
:ю:&A"
 
_output_shapes
:
юм:!B

_output_shapes	
:м:%C!

_output_shapes
:	мd: D

_output_shapes
:d:$E 

_output_shapes

:d2: F

_output_shapes
:2:$G 

_output_shapes

:2: H

_output_shapes
::I

_output_shapes
: 
Їw
Я

G__inference_sequential_1_layer_call_and_return_conditional_losses_76595

inputs&
conv2d_76342:
conv2d_76344:(
conv2d_1_76365:
conv2d_1_76367:(
conv2d_2_76388: 
conv2d_2_76390: (
conv2d_3_76411: @
conv2d_3_76413:@(
conv2d_4_76428:@P
conv2d_4_76430:P
dense_76465:
аю
dense_76467:	ю!
dense_1_76488:
юм
dense_1_76490:	м 
dense_2_76518:	мd
dense_2_76520:d
dense_3_76541:d2
dense_3_76543:2
output_76565:2
output_76567:
identityИвconv2d/StatefulPartitionedCallв conv2d_1/StatefulPartitionedCallв conv2d_2/StatefulPartitionedCallв conv2d_3/StatefulPartitionedCallв conv2d_4/StatefulPartitionedCallвdense/StatefulPartitionedCallв.dense/kernel/Regularizer/Square/ReadVariableOpвdense_1/StatefulPartitionedCallв0dense_1/kernel/Regularizer/Square/ReadVariableOpвdense_2/StatefulPartitionedCallв0dense_2/kernel/Regularizer/Square/ReadVariableOpвdense_3/StatefulPartitionedCallв0dense_3/kernel/Regularizer/Square/ReadVariableOpвoutput/StatefulPartitionedCallс
sequential/PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *N
fIRG
E__inference_sequential_layer_call_and_return_conditional_losses_758282
sequential/PartitionedCallБ
rescaling_1/PartitionedCallPartitionedCall#sequential/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_rescaling_1_layer_call_and_return_conditional_losses_763282
rescaling_1/PartitionedCallн
conv2d/StatefulPartitionedCallStatefulPartitionedCall$rescaling_1/PartitionedCall:output:0conv2d_76342conv2d_76344*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_conv2d_layer_call_and_return_conditional_losses_763412 
conv2d/StatefulPartitionedCallЛ
max_pooling2d/PartitionedCallPartitionedCall'conv2d/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *Q
fLRJ
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_763512
max_pooling2d/PartitionedCall╣
 conv2d_1/StatefulPartitionedCallStatefulPartitionedCall&max_pooling2d/PartitionedCall:output:0conv2d_1_76365conv2d_1_76367*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_1_layer_call_and_return_conditional_losses_763642"
 conv2d_1/StatefulPartitionedCallУ
max_pooling2d_1/PartitionedCallPartitionedCall)conv2d_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_763742!
max_pooling2d_1/PartitionedCall╗
 conv2d_2/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_1/PartitionedCall:output:0conv2d_2_76388conv2d_2_76390*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         	 *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_2_layer_call_and_return_conditional_losses_763872"
 conv2d_2/StatefulPartitionedCallУ
max_pooling2d_2/PartitionedCallPartitionedCall)conv2d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:          * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_763972!
max_pooling2d_2/PartitionedCall╗
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_2/PartitionedCall:output:0conv2d_3_76411conv2d_3_76413*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         @*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_3_layer_call_and_return_conditional_losses_764102"
 conv2d_3/StatefulPartitionedCall╝
 conv2d_4/StatefulPartitionedCallStatefulPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0conv2d_4_76428conv2d_4_76430*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *L
fGRE
C__inference_conv2d_4_layer_call_and_return_conditional_losses_764272"
 conv2d_4/StatefulPartitionedCallУ
max_pooling2d_3/PartitionedCallPartitionedCall)conv2d_4/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         P* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_764372!
max_pooling2d_3/PartitionedCallє
flatten/PartitionedCallPartitionedCall(max_pooling2d_3/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         а* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_764452
flatten/PartitionedCallЭ
dense/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0dense_76465dense_76467*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         ю*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *I
fDRB
@__inference_dense_layer_call_and_return_conditional_losses_764642
dense/StatefulPartitionedCallн
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_76488dense_1_76490*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_1_layer_call_and_return_conditional_losses_764872!
dense_1/StatefulPartitionedCallє
dropout/PartitionedCallPartitionedCall(dense_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:         м* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dropout_layer_call_and_return_conditional_losses_764982
dropout/PartitionedCallж
dense_2/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0dense_2_76518dense_2_76520*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         d*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_765172!
dense_2/StatefulPartitionedCallо
dense_3/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0dense_3_76541dense_3_76543*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_765402!
dense_3/StatefulPartitionedCall°
dropout_1/PartitionedCallPartitionedCall(dense_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         2* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *M
fHRF
D__inference_dropout_1_layer_call_and_return_conditional_losses_765512
dropout_1/PartitionedCallг
output/StatefulPartitionedCallStatefulPartitionedCall"dropout_1/PartitionedCall:output:0output_76565output_76567*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_output_layer_call_and_return_conditional_losses_765642 
output/StatefulPartitionedCallо
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_76465* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul┤
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_76488* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul│
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_76518*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul▓
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_76541*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mulВ
IdentityIdentity'output/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identityь
NoOpNoOp^conv2d/StatefulPartitionedCall!^conv2d_1/StatefulPartitionedCall!^conv2d_2/StatefulPartitionedCall!^conv2d_3/StatefulPartitionedCall!^conv2d_4/StatefulPartitionedCall^dense/StatefulPartitionedCall/^dense/kernel/Regularizer/Square/ReadVariableOp ^dense_1/StatefulPartitionedCall1^dense_1/kernel/Regularizer/Square/ReadVariableOp ^dense_2/StatefulPartitionedCall1^dense_2/kernel/Regularizer/Square/ReadVariableOp ^dense_3/StatefulPartitionedCall1^dense_3/kernel/Regularizer/Square/ReadVariableOp^output/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*V
_input_shapesE
C:         &: : : : : : : : : : : : : : : : : : : : 2@
conv2d/StatefulPartitionedCallconv2d/StatefulPartitionedCall2D
 conv2d_1/StatefulPartitionedCall conv2d_1/StatefulPartitionedCall2D
 conv2d_2/StatefulPartitionedCall conv2d_2/StatefulPartitionedCall2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2D
 conv2d_4/StatefulPartitionedCall conv2d_4/StatefulPartitionedCall2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2@
output/StatefulPartitionedCalloutput/StatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
щ
ж
B__inference_dense_3_layer_call_and_return_conditional_losses_76540

inputs0
matmul_readvariableop_resource:d2-
biasadd_readvariableop_resource:2
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpв0dense_3/kernel/Regularizer/Square/ReadVariableOpН
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:d2*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:         22
Relu├
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mulm
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:         22

Identity▓
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         d: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:         d
 
_user_specified_nameinputs
ч
·
A__inference_conv2d_layer_call_and_return_conditional_losses_76341

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &*
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         &2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         &2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         &: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
Ь▓
Ї+
!__inference__traced_restore_79503
file_prefix8
assignvariableop_conv2d_kernel:,
assignvariableop_1_conv2d_bias:<
"assignvariableop_2_conv2d_1_kernel:.
 assignvariableop_3_conv2d_1_bias:<
"assignvariableop_4_conv2d_2_kernel: .
 assignvariableop_5_conv2d_2_bias: <
"assignvariableop_6_conv2d_3_kernel: @.
 assignvariableop_7_conv2d_3_bias:@<
"assignvariableop_8_conv2d_4_kernel:@P.
 assignvariableop_9_conv2d_4_bias:P4
 assignvariableop_10_dense_kernel:
аю-
assignvariableop_11_dense_bias:	ю6
"assignvariableop_12_dense_1_kernel:
юм/
 assignvariableop_13_dense_1_bias:	м5
"assignvariableop_14_dense_2_kernel:	мd.
 assignvariableop_15_dense_2_bias:d4
"assignvariableop_16_dense_3_kernel:d2.
 assignvariableop_17_dense_3_bias:23
!assignvariableop_18_output_kernel:2-
assignvariableop_19_output_bias:'
assignvariableop_20_adam_iter:	 )
assignvariableop_21_adam_beta_1: )
assignvariableop_22_adam_beta_2: (
assignvariableop_23_adam_decay: 0
&assignvariableop_24_adam_learning_rate: *
assignvariableop_25_variable:	,
assignvariableop_26_variable_1:	,
assignvariableop_27_variable_2:	#
assignvariableop_28_total: #
assignvariableop_29_count: %
assignvariableop_30_total_1: %
assignvariableop_31_count_1: B
(assignvariableop_32_adam_conv2d_kernel_m:4
&assignvariableop_33_adam_conv2d_bias_m:D
*assignvariableop_34_adam_conv2d_1_kernel_m:6
(assignvariableop_35_adam_conv2d_1_bias_m:D
*assignvariableop_36_adam_conv2d_2_kernel_m: 6
(assignvariableop_37_adam_conv2d_2_bias_m: D
*assignvariableop_38_adam_conv2d_3_kernel_m: @6
(assignvariableop_39_adam_conv2d_3_bias_m:@D
*assignvariableop_40_adam_conv2d_4_kernel_m:@P6
(assignvariableop_41_adam_conv2d_4_bias_m:P;
'assignvariableop_42_adam_dense_kernel_m:
аю4
%assignvariableop_43_adam_dense_bias_m:	ю=
)assignvariableop_44_adam_dense_1_kernel_m:
юм6
'assignvariableop_45_adam_dense_1_bias_m:	м<
)assignvariableop_46_adam_dense_2_kernel_m:	мd5
'assignvariableop_47_adam_dense_2_bias_m:d;
)assignvariableop_48_adam_dense_3_kernel_m:d25
'assignvariableop_49_adam_dense_3_bias_m:2:
(assignvariableop_50_adam_output_kernel_m:24
&assignvariableop_51_adam_output_bias_m:B
(assignvariableop_52_adam_conv2d_kernel_v:4
&assignvariableop_53_adam_conv2d_bias_v:D
*assignvariableop_54_adam_conv2d_1_kernel_v:6
(assignvariableop_55_adam_conv2d_1_bias_v:D
*assignvariableop_56_adam_conv2d_2_kernel_v: 6
(assignvariableop_57_adam_conv2d_2_bias_v: D
*assignvariableop_58_adam_conv2d_3_kernel_v: @6
(assignvariableop_59_adam_conv2d_3_bias_v:@D
*assignvariableop_60_adam_conv2d_4_kernel_v:@P6
(assignvariableop_61_adam_conv2d_4_bias_v:P;
'assignvariableop_62_adam_dense_kernel_v:
аю4
%assignvariableop_63_adam_dense_bias_v:	ю=
)assignvariableop_64_adam_dense_1_kernel_v:
юм6
'assignvariableop_65_adam_dense_1_bias_v:	м<
)assignvariableop_66_adam_dense_2_kernel_v:	мd5
'assignvariableop_67_adam_dense_2_bias_v:d;
)assignvariableop_68_adam_dense_3_kernel_v:d25
'assignvariableop_69_adam_dense_3_bias_v:2:
(assignvariableop_70_adam_output_kernel_v:24
&assignvariableop_71_adam_output_bias_v:
identity_73ИвAssignVariableOpвAssignVariableOp_1вAssignVariableOp_10вAssignVariableOp_11вAssignVariableOp_12вAssignVariableOp_13вAssignVariableOp_14вAssignVariableOp_15вAssignVariableOp_16вAssignVariableOp_17вAssignVariableOp_18вAssignVariableOp_19вAssignVariableOp_2вAssignVariableOp_20вAssignVariableOp_21вAssignVariableOp_22вAssignVariableOp_23вAssignVariableOp_24вAssignVariableOp_25вAssignVariableOp_26вAssignVariableOp_27вAssignVariableOp_28вAssignVariableOp_29вAssignVariableOp_3вAssignVariableOp_30вAssignVariableOp_31вAssignVariableOp_32вAssignVariableOp_33вAssignVariableOp_34вAssignVariableOp_35вAssignVariableOp_36вAssignVariableOp_37вAssignVariableOp_38вAssignVariableOp_39вAssignVariableOp_4вAssignVariableOp_40вAssignVariableOp_41вAssignVariableOp_42вAssignVariableOp_43вAssignVariableOp_44вAssignVariableOp_45вAssignVariableOp_46вAssignVariableOp_47вAssignVariableOp_48вAssignVariableOp_49вAssignVariableOp_5вAssignVariableOp_50вAssignVariableOp_51вAssignVariableOp_52вAssignVariableOp_53вAssignVariableOp_54вAssignVariableOp_55вAssignVariableOp_56вAssignVariableOp_57вAssignVariableOp_58вAssignVariableOp_59вAssignVariableOp_6вAssignVariableOp_60вAssignVariableOp_61вAssignVariableOp_62вAssignVariableOp_63вAssignVariableOp_64вAssignVariableOp_65вAssignVariableOp_66вAssignVariableOp_67вAssignVariableOp_68вAssignVariableOp_69вAssignVariableOp_7вAssignVariableOp_70вAssignVariableOp_71вAssignVariableOp_8вAssignVariableOp_9▄(
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:I*
dtype0*ш'
value▐'B█'IB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-9/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-9/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB:layer-0/layer-0/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEB:layer-0/layer-1/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEB:layer-0/layer-2/_rng/_state_var/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-9/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-9/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-4/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-4/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-5/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-5/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-6/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-6/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-7/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-7/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-8/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-8/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-9/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-9/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_namesг
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:I*
dtype0*з
valueЭBЪIB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slicesЫ
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*║
_output_shapesз
д:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*W
dtypesM
K2I				2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

IdentityЭ
AssignVariableOpAssignVariableOpassignvariableop_conv2d_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1г
AssignVariableOp_1AssignVariableOpassignvariableop_1_conv2d_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2з
AssignVariableOp_2AssignVariableOp"assignvariableop_2_conv2d_1_kernelIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3е
AssignVariableOp_3AssignVariableOp assignvariableop_3_conv2d_1_biasIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4з
AssignVariableOp_4AssignVariableOp"assignvariableop_4_conv2d_2_kernelIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5е
AssignVariableOp_5AssignVariableOp assignvariableop_5_conv2d_2_biasIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6з
AssignVariableOp_6AssignVariableOp"assignvariableop_6_conv2d_3_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7е
AssignVariableOp_7AssignVariableOp assignvariableop_7_conv2d_3_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8з
AssignVariableOp_8AssignVariableOp"assignvariableop_8_conv2d_4_kernelIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9е
AssignVariableOp_9AssignVariableOp assignvariableop_9_conv2d_4_biasIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10и
AssignVariableOp_10AssignVariableOp assignvariableop_10_dense_kernelIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11ж
AssignVariableOp_11AssignVariableOpassignvariableop_11_dense_biasIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:2
Identity_12к
AssignVariableOp_12AssignVariableOp"assignvariableop_12_dense_1_kernelIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13и
AssignVariableOp_13AssignVariableOp assignvariableop_13_dense_1_biasIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14к
AssignVariableOp_14AssignVariableOp"assignvariableop_14_dense_2_kernelIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15и
AssignVariableOp_15AssignVariableOp assignvariableop_15_dense_2_biasIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16к
AssignVariableOp_16AssignVariableOp"assignvariableop_16_dense_3_kernelIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17и
AssignVariableOp_17AssignVariableOp assignvariableop_17_dense_3_biasIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18й
AssignVariableOp_18AssignVariableOp!assignvariableop_18_output_kernelIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19з
AssignVariableOp_19AssignVariableOpassignvariableop_19_output_biasIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_20е
AssignVariableOp_20AssignVariableOpassignvariableop_20_adam_iterIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21з
AssignVariableOp_21AssignVariableOpassignvariableop_21_adam_beta_1Identity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22з
AssignVariableOp_22AssignVariableOpassignvariableop_22_adam_beta_2Identity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23ж
AssignVariableOp_23AssignVariableOpassignvariableop_23_adam_decayIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24о
AssignVariableOp_24AssignVariableOp&assignvariableop_24_adam_learning_rateIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_25д
AssignVariableOp_25AssignVariableOpassignvariableop_25_variableIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_26ж
AssignVariableOp_26AssignVariableOpassignvariableop_26_variable_1Identity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0	*
_output_shapes
:2
Identity_27ж
AssignVariableOp_27AssignVariableOpassignvariableop_27_variable_2Identity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28б
AssignVariableOp_28AssignVariableOpassignvariableop_28_totalIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29б
AssignVariableOp_29AssignVariableOpassignvariableop_29_countIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30г
AssignVariableOp_30AssignVariableOpassignvariableop_30_total_1Identity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31г
AssignVariableOp_31AssignVariableOpassignvariableop_31_count_1Identity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32░
AssignVariableOp_32AssignVariableOp(assignvariableop_32_adam_conv2d_kernel_mIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33о
AssignVariableOp_33AssignVariableOp&assignvariableop_33_adam_conv2d_bias_mIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_33n
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:2
Identity_34▓
AssignVariableOp_34AssignVariableOp*assignvariableop_34_adam_conv2d_1_kernel_mIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_34n
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:2
Identity_35░
AssignVariableOp_35AssignVariableOp(assignvariableop_35_adam_conv2d_1_bias_mIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_35n
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:2
Identity_36▓
AssignVariableOp_36AssignVariableOp*assignvariableop_36_adam_conv2d_2_kernel_mIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_36n
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:2
Identity_37░
AssignVariableOp_37AssignVariableOp(assignvariableop_37_adam_conv2d_2_bias_mIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_37n
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:2
Identity_38▓
AssignVariableOp_38AssignVariableOp*assignvariableop_38_adam_conv2d_3_kernel_mIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_38n
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:2
Identity_39░
AssignVariableOp_39AssignVariableOp(assignvariableop_39_adam_conv2d_3_bias_mIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_39n
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:2
Identity_40▓
AssignVariableOp_40AssignVariableOp*assignvariableop_40_adam_conv2d_4_kernel_mIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_40n
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:2
Identity_41░
AssignVariableOp_41AssignVariableOp(assignvariableop_41_adam_conv2d_4_bias_mIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_41n
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:2
Identity_42п
AssignVariableOp_42AssignVariableOp'assignvariableop_42_adam_dense_kernel_mIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_42n
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:2
Identity_43н
AssignVariableOp_43AssignVariableOp%assignvariableop_43_adam_dense_bias_mIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_43n
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:2
Identity_44▒
AssignVariableOp_44AssignVariableOp)assignvariableop_44_adam_dense_1_kernel_mIdentity_44:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_44n
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:2
Identity_45п
AssignVariableOp_45AssignVariableOp'assignvariableop_45_adam_dense_1_bias_mIdentity_45:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_45n
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:2
Identity_46▒
AssignVariableOp_46AssignVariableOp)assignvariableop_46_adam_dense_2_kernel_mIdentity_46:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_46n
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:2
Identity_47п
AssignVariableOp_47AssignVariableOp'assignvariableop_47_adam_dense_2_bias_mIdentity_47:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_47n
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:2
Identity_48▒
AssignVariableOp_48AssignVariableOp)assignvariableop_48_adam_dense_3_kernel_mIdentity_48:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_48n
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:2
Identity_49п
AssignVariableOp_49AssignVariableOp'assignvariableop_49_adam_dense_3_bias_mIdentity_49:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_49n
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:2
Identity_50░
AssignVariableOp_50AssignVariableOp(assignvariableop_50_adam_output_kernel_mIdentity_50:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_50n
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:2
Identity_51о
AssignVariableOp_51AssignVariableOp&assignvariableop_51_adam_output_bias_mIdentity_51:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_51n
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:2
Identity_52░
AssignVariableOp_52AssignVariableOp(assignvariableop_52_adam_conv2d_kernel_vIdentity_52:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_52n
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:2
Identity_53о
AssignVariableOp_53AssignVariableOp&assignvariableop_53_adam_conv2d_bias_vIdentity_53:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_53n
Identity_54IdentityRestoreV2:tensors:54"/device:CPU:0*
T0*
_output_shapes
:2
Identity_54▓
AssignVariableOp_54AssignVariableOp*assignvariableop_54_adam_conv2d_1_kernel_vIdentity_54:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_54n
Identity_55IdentityRestoreV2:tensors:55"/device:CPU:0*
T0*
_output_shapes
:2
Identity_55░
AssignVariableOp_55AssignVariableOp(assignvariableop_55_adam_conv2d_1_bias_vIdentity_55:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_55n
Identity_56IdentityRestoreV2:tensors:56"/device:CPU:0*
T0*
_output_shapes
:2
Identity_56▓
AssignVariableOp_56AssignVariableOp*assignvariableop_56_adam_conv2d_2_kernel_vIdentity_56:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_56n
Identity_57IdentityRestoreV2:tensors:57"/device:CPU:0*
T0*
_output_shapes
:2
Identity_57░
AssignVariableOp_57AssignVariableOp(assignvariableop_57_adam_conv2d_2_bias_vIdentity_57:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_57n
Identity_58IdentityRestoreV2:tensors:58"/device:CPU:0*
T0*
_output_shapes
:2
Identity_58▓
AssignVariableOp_58AssignVariableOp*assignvariableop_58_adam_conv2d_3_kernel_vIdentity_58:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_58n
Identity_59IdentityRestoreV2:tensors:59"/device:CPU:0*
T0*
_output_shapes
:2
Identity_59░
AssignVariableOp_59AssignVariableOp(assignvariableop_59_adam_conv2d_3_bias_vIdentity_59:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_59n
Identity_60IdentityRestoreV2:tensors:60"/device:CPU:0*
T0*
_output_shapes
:2
Identity_60▓
AssignVariableOp_60AssignVariableOp*assignvariableop_60_adam_conv2d_4_kernel_vIdentity_60:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_60n
Identity_61IdentityRestoreV2:tensors:61"/device:CPU:0*
T0*
_output_shapes
:2
Identity_61░
AssignVariableOp_61AssignVariableOp(assignvariableop_61_adam_conv2d_4_bias_vIdentity_61:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_61n
Identity_62IdentityRestoreV2:tensors:62"/device:CPU:0*
T0*
_output_shapes
:2
Identity_62п
AssignVariableOp_62AssignVariableOp'assignvariableop_62_adam_dense_kernel_vIdentity_62:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_62n
Identity_63IdentityRestoreV2:tensors:63"/device:CPU:0*
T0*
_output_shapes
:2
Identity_63н
AssignVariableOp_63AssignVariableOp%assignvariableop_63_adam_dense_bias_vIdentity_63:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_63n
Identity_64IdentityRestoreV2:tensors:64"/device:CPU:0*
T0*
_output_shapes
:2
Identity_64▒
AssignVariableOp_64AssignVariableOp)assignvariableop_64_adam_dense_1_kernel_vIdentity_64:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_64n
Identity_65IdentityRestoreV2:tensors:65"/device:CPU:0*
T0*
_output_shapes
:2
Identity_65п
AssignVariableOp_65AssignVariableOp'assignvariableop_65_adam_dense_1_bias_vIdentity_65:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_65n
Identity_66IdentityRestoreV2:tensors:66"/device:CPU:0*
T0*
_output_shapes
:2
Identity_66▒
AssignVariableOp_66AssignVariableOp)assignvariableop_66_adam_dense_2_kernel_vIdentity_66:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_66n
Identity_67IdentityRestoreV2:tensors:67"/device:CPU:0*
T0*
_output_shapes
:2
Identity_67п
AssignVariableOp_67AssignVariableOp'assignvariableop_67_adam_dense_2_bias_vIdentity_67:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_67n
Identity_68IdentityRestoreV2:tensors:68"/device:CPU:0*
T0*
_output_shapes
:2
Identity_68▒
AssignVariableOp_68AssignVariableOp)assignvariableop_68_adam_dense_3_kernel_vIdentity_68:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_68n
Identity_69IdentityRestoreV2:tensors:69"/device:CPU:0*
T0*
_output_shapes
:2
Identity_69п
AssignVariableOp_69AssignVariableOp'assignvariableop_69_adam_dense_3_bias_vIdentity_69:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_69n
Identity_70IdentityRestoreV2:tensors:70"/device:CPU:0*
T0*
_output_shapes
:2
Identity_70░
AssignVariableOp_70AssignVariableOp(assignvariableop_70_adam_output_kernel_vIdentity_70:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_70n
Identity_71IdentityRestoreV2:tensors:71"/device:CPU:0*
T0*
_output_shapes
:2
Identity_71о
AssignVariableOp_71AssignVariableOp&assignvariableop_71_adam_output_bias_vIdentity_71:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_719
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOpО
Identity_72Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_62^AssignVariableOp_63^AssignVariableOp_64^AssignVariableOp_65^AssignVariableOp_66^AssignVariableOp_67^AssignVariableOp_68^AssignVariableOp_69^AssignVariableOp_7^AssignVariableOp_70^AssignVariableOp_71^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_72f
Identity_73IdentityIdentity_72:output:0^NoOp_1*
T0*
_output_shapes
: 2
Identity_73Ў
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_62^AssignVariableOp_63^AssignVariableOp_64^AssignVariableOp_65^AssignVariableOp_66^AssignVariableOp_67^AssignVariableOp_68^AssignVariableOp_69^AssignVariableOp_7^AssignVariableOp_70^AssignVariableOp_71^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 2
NoOp_1"#
identity_73Identity_73:output:0*з
_input_shapesХ
Т: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2$
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
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522*
AssignVariableOp_53AssignVariableOp_532*
AssignVariableOp_54AssignVariableOp_542*
AssignVariableOp_55AssignVariableOp_552*
AssignVariableOp_56AssignVariableOp_562*
AssignVariableOp_57AssignVariableOp_572*
AssignVariableOp_58AssignVariableOp_582*
AssignVariableOp_59AssignVariableOp_592(
AssignVariableOp_6AssignVariableOp_62*
AssignVariableOp_60AssignVariableOp_602*
AssignVariableOp_61AssignVariableOp_612*
AssignVariableOp_62AssignVariableOp_622*
AssignVariableOp_63AssignVariableOp_632*
AssignVariableOp_64AssignVariableOp_642*
AssignVariableOp_65AssignVariableOp_652*
AssignVariableOp_66AssignVariableOp_662*
AssignVariableOp_67AssignVariableOp_672*
AssignVariableOp_68AssignVariableOp_682*
AssignVariableOp_69AssignVariableOp_692(
AssignVariableOp_7AssignVariableOp_72*
AssignVariableOp_70AssignVariableOp_702*
AssignVariableOp_71AssignVariableOp_712(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
Ъ
b
F__inference_random_flip_layer_call_and_return_conditional_losses_75813

inputs
identityb
IdentityIdentityinputs*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
и
d
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_78279

inputs
identityн
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4                                    *
ksize
*
paddingVALID*
strides
2	
MaxPoolЗ
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
тх
▐
G__inference_sequential_1_layer_call_and_return_conditional_losses_77925

inputsV
Hsequential_random_flip_stateful_uniform_full_int_rngreadandskip_resource:	Q
Csequential_random_rotation_stateful_uniform_rngreadandskip_resource:	M
?sequential_random_zoom_stateful_uniform_rngreadandskip_resource:	?
%conv2d_conv2d_readvariableop_resource:4
&conv2d_biasadd_readvariableop_resource:A
'conv2d_1_conv2d_readvariableop_resource:6
(conv2d_1_biasadd_readvariableop_resource:A
'conv2d_2_conv2d_readvariableop_resource: 6
(conv2d_2_biasadd_readvariableop_resource: A
'conv2d_3_conv2d_readvariableop_resource: @6
(conv2d_3_biasadd_readvariableop_resource:@A
'conv2d_4_conv2d_readvariableop_resource:@P6
(conv2d_4_biasadd_readvariableop_resource:P8
$dense_matmul_readvariableop_resource:
аю4
%dense_biasadd_readvariableop_resource:	ю:
&dense_1_matmul_readvariableop_resource:
юм6
'dense_1_biasadd_readvariableop_resource:	м9
&dense_2_matmul_readvariableop_resource:	мd5
'dense_2_biasadd_readvariableop_resource:d8
&dense_3_matmul_readvariableop_resource:d25
'dense_3_biasadd_readvariableop_resource:27
%output_matmul_readvariableop_resource:24
&output_biasadd_readvariableop_resource:
identityИвconv2d/BiasAdd/ReadVariableOpвconv2d/Conv2D/ReadVariableOpвconv2d_1/BiasAdd/ReadVariableOpвconv2d_1/Conv2D/ReadVariableOpвconv2d_2/BiasAdd/ReadVariableOpвconv2d_2/Conv2D/ReadVariableOpвconv2d_3/BiasAdd/ReadVariableOpвconv2d_3/Conv2D/ReadVariableOpвconv2d_4/BiasAdd/ReadVariableOpвconv2d_4/Conv2D/ReadVariableOpвdense/BiasAdd/ReadVariableOpвdense/MatMul/ReadVariableOpв.dense/kernel/Regularizer/Square/ReadVariableOpвdense_1/BiasAdd/ReadVariableOpвdense_1/MatMul/ReadVariableOpв0dense_1/kernel/Regularizer/Square/ReadVariableOpвdense_2/BiasAdd/ReadVariableOpвdense_2/MatMul/ReadVariableOpв0dense_2/kernel/Regularizer/Square/ReadVariableOpвdense_3/BiasAdd/ReadVariableOpвdense_3/MatMul/ReadVariableOpв0dense_3/kernel/Regularizer/Square/ReadVariableOpвoutput/BiasAdd/ReadVariableOpвoutput/MatMul/ReadVariableOpв?sequential/random_flip/stateful_uniform_full_int/RngReadAndSkipвfsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgвmsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterв:sequential/random_rotation/stateful_uniform/RngReadAndSkipв6sequential/random_zoom/stateful_uniform/RngReadAndSkip║
6sequential/random_flip/stateful_uniform_full_int/shapeConst*
_output_shapes
:*
dtype0*
valueB:28
6sequential/random_flip/stateful_uniform_full_int/shape║
6sequential/random_flip/stateful_uniform_full_int/ConstConst*
_output_shapes
:*
dtype0*
valueB: 28
6sequential/random_flip/stateful_uniform_full_int/ConstЩ
5sequential/random_flip/stateful_uniform_full_int/ProdProd?sequential/random_flip/stateful_uniform_full_int/shape:output:0?sequential/random_flip/stateful_uniform_full_int/Const:output:0*
T0*
_output_shapes
: 27
5sequential/random_flip/stateful_uniform_full_int/Prod┤
7sequential/random_flip/stateful_uniform_full_int/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :29
7sequential/random_flip/stateful_uniform_full_int/Cast/xъ
7sequential/random_flip/stateful_uniform_full_int/Cast_1Cast>sequential/random_flip/stateful_uniform_full_int/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 29
7sequential/random_flip/stateful_uniform_full_int/Cast_1∙
?sequential/random_flip/stateful_uniform_full_int/RngReadAndSkipRngReadAndSkipHsequential_random_flip_stateful_uniform_full_int_rngreadandskip_resource@sequential/random_flip/stateful_uniform_full_int/Cast/x:output:0;sequential/random_flip/stateful_uniform_full_int/Cast_1:y:0*
_output_shapes
:2A
?sequential/random_flip/stateful_uniform_full_int/RngReadAndSkip╓
Dsequential/random_flip/stateful_uniform_full_int/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2F
Dsequential/random_flip/stateful_uniform_full_int/strided_slice/stack┌
Fsequential/random_flip/stateful_uniform_full_int/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2H
Fsequential/random_flip/stateful_uniform_full_int/strided_slice/stack_1┌
Fsequential/random_flip/stateful_uniform_full_int/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2H
Fsequential/random_flip/stateful_uniform_full_int/strided_slice/stack_2О
>sequential/random_flip/stateful_uniform_full_int/strided_sliceStridedSliceGsequential/random_flip/stateful_uniform_full_int/RngReadAndSkip:value:0Msequential/random_flip/stateful_uniform_full_int/strided_slice/stack:output:0Osequential/random_flip/stateful_uniform_full_int/strided_slice/stack_1:output:0Osequential/random_flip/stateful_uniform_full_int/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2@
>sequential/random_flip/stateful_uniform_full_int/strided_slice∙
8sequential/random_flip/stateful_uniform_full_int/BitcastBitcastGsequential/random_flip/stateful_uniform_full_int/strided_slice:output:0*
T0	*
_output_shapes
:*

type02:
8sequential/random_flip/stateful_uniform_full_int/Bitcast┌
Fsequential/random_flip/stateful_uniform_full_int/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2H
Fsequential/random_flip/stateful_uniform_full_int/strided_slice_1/stack▐
Hsequential/random_flip/stateful_uniform_full_int/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2J
Hsequential/random_flip/stateful_uniform_full_int/strided_slice_1/stack_1▐
Hsequential/random_flip/stateful_uniform_full_int/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2J
Hsequential/random_flip/stateful_uniform_full_int/strided_slice_1/stack_2Ж
@sequential/random_flip/stateful_uniform_full_int/strided_slice_1StridedSliceGsequential/random_flip/stateful_uniform_full_int/RngReadAndSkip:value:0Osequential/random_flip/stateful_uniform_full_int/strided_slice_1/stack:output:0Qsequential/random_flip/stateful_uniform_full_int/strided_slice_1/stack_1:output:0Qsequential/random_flip/stateful_uniform_full_int/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2B
@sequential/random_flip/stateful_uniform_full_int/strided_slice_1 
:sequential/random_flip/stateful_uniform_full_int/Bitcast_1BitcastIsequential/random_flip/stateful_uniform_full_int/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02<
:sequential/random_flip/stateful_uniform_full_int/Bitcast_1о
4sequential/random_flip/stateful_uniform_full_int/algConst*
_output_shapes
: *
dtype0*
value	B :26
4sequential/random_flip/stateful_uniform_full_int/alg╕
0sequential/random_flip/stateful_uniform_full_intStatelessRandomUniformFullIntV2?sequential/random_flip/stateful_uniform_full_int/shape:output:0Csequential/random_flip/stateful_uniform_full_int/Bitcast_1:output:0Asequential/random_flip/stateful_uniform_full_int/Bitcast:output:0=sequential/random_flip/stateful_uniform_full_int/alg:output:0*
_output_shapes
:*
dtype0	22
0sequential/random_flip/stateful_uniform_full_intР
!sequential/random_flip/zeros_likeConst*
_output_shapes
:*
dtype0	*
valueB	R 2#
!sequential/random_flip/zeros_like▌
sequential/random_flip/stackPack9sequential/random_flip/stateful_uniform_full_int:output:0*sequential/random_flip/zeros_like:output:0*
N*
T0	*
_output_shapes

:2
sequential/random_flip/stackй
*sequential/random_flip/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB"        2,
*sequential/random_flip/strided_slice/stackн
,sequential/random_flip/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB"       2.
,sequential/random_flip/strided_slice/stack_1н
,sequential/random_flip/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2.
,sequential/random_flip/strided_slice/stack_2Т
$sequential/random_flip/strided_sliceStridedSlice%sequential/random_flip/stack:output:03sequential/random_flip/strided_slice/stack:output:05sequential/random_flip/strided_slice/stack_1:output:05sequential/random_flip/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask*
end_mask*
shrink_axis_mask2&
$sequential/random_flip/strided_sliceБ
Jsequential/random_flip/stateless_random_flip_left_right/control_dependencyIdentityinputs*
T0*
_class
loc:@inputs*/
_output_shapes
:         &2L
Jsequential/random_flip/stateless_random_flip_left_right/control_dependencyБ
=sequential/random_flip/stateless_random_flip_left_right/ShapeShapeSsequential/random_flip/stateless_random_flip_left_right/control_dependency:output:0*
T0*
_output_shapes
:2?
=sequential/random_flip/stateless_random_flip_left_right/Shapeф
Ksequential/random_flip/stateless_random_flip_left_right/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2M
Ksequential/random_flip/stateless_random_flip_left_right/strided_slice/stackш
Msequential/random_flip/stateless_random_flip_left_right/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2O
Msequential/random_flip/stateless_random_flip_left_right/strided_slice/stack_1ш
Msequential/random_flip/stateless_random_flip_left_right/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2O
Msequential/random_flip/stateless_random_flip_left_right/strided_slice/stack_2▓
Esequential/random_flip/stateless_random_flip_left_right/strided_sliceStridedSliceFsequential/random_flip/stateless_random_flip_left_right/Shape:output:0Tsequential/random_flip/stateless_random_flip_left_right/strided_slice/stack:output:0Vsequential/random_flip/stateless_random_flip_left_right/strided_slice/stack_1:output:0Vsequential/random_flip/stateless_random_flip_left_right/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2G
Esequential/random_flip/stateless_random_flip_left_right/strided_slice╢
Vsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/shapePackNsequential/random_flip/stateless_random_flip_left_right/strided_slice:output:0*
N*
T0*
_output_shapes
:2X
Vsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/shapeё
Tsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *    2V
Tsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/minё
Tsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2V
Tsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/max╧
msequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterStatelessRandomGetKeyCounter-sequential/random_flip/strided_slice:output:0* 
_output_shapes
::2o
msequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounterё
fsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgStatelessRandomGetAlgn^sequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter*
_output_shapes
: 2h
fsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg╘
isequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2StatelessRandomUniformV2_sequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/shape:output:0ssequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:key:0wsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter:counter:0lsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg:alg:0*#
_output_shapes
:         2k
isequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2Т
Tsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/subSub]sequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/max:output:0]sequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/min:output:0*
T0*
_output_shapes
: 2V
Tsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/subп
Tsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/mulMulrsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomUniformV2:output:0Xsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/sub:z:0*
T0*#
_output_shapes
:         2V
Tsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/mulФ
Psequential/random_flip/stateless_random_flip_left_right/stateless_random_uniformAddV2Xsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/mul:z:0]sequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/min:output:0*
T0*#
_output_shapes
:         2R
Psequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform╘
Gsequential/random_flip/stateless_random_flip_left_right/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2I
Gsequential/random_flip/stateless_random_flip_left_right/Reshape/shape/1╘
Gsequential/random_flip/stateless_random_flip_left_right/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2I
Gsequential/random_flip/stateless_random_flip_left_right/Reshape/shape/2╘
Gsequential/random_flip/stateless_random_flip_left_right/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2I
Gsequential/random_flip/stateless_random_flip_left_right/Reshape/shape/3К
Esequential/random_flip/stateless_random_flip_left_right/Reshape/shapePackNsequential/random_flip/stateless_random_flip_left_right/strided_slice:output:0Psequential/random_flip/stateless_random_flip_left_right/Reshape/shape/1:output:0Psequential/random_flip/stateless_random_flip_left_right/Reshape/shape/2:output:0Psequential/random_flip/stateless_random_flip_left_right/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2G
Esequential/random_flip/stateless_random_flip_left_right/Reshape/shapeэ
?sequential/random_flip/stateless_random_flip_left_right/ReshapeReshapeTsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform:z:0Nsequential/random_flip/stateless_random_flip_left_right/Reshape/shape:output:0*
T0*/
_output_shapes
:         2A
?sequential/random_flip/stateless_random_flip_left_right/ReshapeЛ
=sequential/random_flip/stateless_random_flip_left_right/RoundRoundHsequential/random_flip/stateless_random_flip_left_right/Reshape:output:0*
T0*/
_output_shapes
:         2?
=sequential/random_flip/stateless_random_flip_left_right/Round┌
Fsequential/random_flip/stateless_random_flip_left_right/ReverseV2/axisConst*
_output_shapes
:*
dtype0*
valueB:2H
Fsequential/random_flip/stateless_random_flip_left_right/ReverseV2/axisє
Asequential/random_flip/stateless_random_flip_left_right/ReverseV2	ReverseV2Ssequential/random_flip/stateless_random_flip_left_right/control_dependency:output:0Osequential/random_flip/stateless_random_flip_left_right/ReverseV2/axis:output:0*
T0*/
_output_shapes
:         &2C
Asequential/random_flip/stateless_random_flip_left_right/ReverseV2╩
;sequential/random_flip/stateless_random_flip_left_right/mulMulAsequential/random_flip/stateless_random_flip_left_right/Round:y:0Jsequential/random_flip/stateless_random_flip_left_right/ReverseV2:output:0*
T0*/
_output_shapes
:         &2=
;sequential/random_flip/stateless_random_flip_left_right/mul├
=sequential/random_flip/stateless_random_flip_left_right/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2?
=sequential/random_flip/stateless_random_flip_left_right/sub/x╞
;sequential/random_flip/stateless_random_flip_left_right/subSubFsequential/random_flip/stateless_random_flip_left_right/sub/x:output:0Asequential/random_flip/stateless_random_flip_left_right/Round:y:0*
T0*/
_output_shapes
:         2=
;sequential/random_flip/stateless_random_flip_left_right/sub╒
=sequential/random_flip/stateless_random_flip_left_right/mul_1Mul?sequential/random_flip/stateless_random_flip_left_right/sub:z:0Ssequential/random_flip/stateless_random_flip_left_right/control_dependency:output:0*
T0*/
_output_shapes
:         &2?
=sequential/random_flip/stateless_random_flip_left_right/mul_1┴
;sequential/random_flip/stateless_random_flip_left_right/addAddV2?sequential/random_flip/stateless_random_flip_left_right/mul:z:0Asequential/random_flip/stateless_random_flip_left_right/mul_1:z:0*
T0*/
_output_shapes
:         &2=
;sequential/random_flip/stateless_random_flip_left_right/add│
 sequential/random_rotation/ShapeShape?sequential/random_flip/stateless_random_flip_left_right/add:z:0*
T0*
_output_shapes
:2"
 sequential/random_rotation/Shapeк
.sequential/random_rotation/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 20
.sequential/random_rotation/strided_slice/stackо
0sequential/random_rotation/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:22
0sequential/random_rotation/strided_slice/stack_1о
0sequential/random_rotation/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:22
0sequential/random_rotation/strided_slice/stack_2Д
(sequential/random_rotation/strided_sliceStridedSlice)sequential/random_rotation/Shape:output:07sequential/random_rotation/strided_slice/stack:output:09sequential/random_rotation/strided_slice/stack_1:output:09sequential/random_rotation/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2*
(sequential/random_rotation/strided_slice╖
0sequential/random_rotation/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:
¤        22
0sequential/random_rotation/strided_slice_1/stack╗
2sequential/random_rotation/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
■        24
2sequential/random_rotation/strided_slice_1/stack_1▓
2sequential/random_rotation/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:24
2sequential/random_rotation/strided_slice_1/stack_2О
*sequential/random_rotation/strided_slice_1StridedSlice)sequential/random_rotation/Shape:output:09sequential/random_rotation/strided_slice_1/stack:output:0;sequential/random_rotation/strided_slice_1/stack_1:output:0;sequential/random_rotation/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2,
*sequential/random_rotation/strided_slice_1п
sequential/random_rotation/CastCast3sequential/random_rotation/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2!
sequential/random_rotation/Cast╖
0sequential/random_rotation/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        22
0sequential/random_rotation/strided_slice_2/stack╗
2sequential/random_rotation/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         24
2sequential/random_rotation/strided_slice_2/stack_1▓
2sequential/random_rotation/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:24
2sequential/random_rotation/strided_slice_2/stack_2О
*sequential/random_rotation/strided_slice_2StridedSlice)sequential/random_rotation/Shape:output:09sequential/random_rotation/strided_slice_2/stack:output:0;sequential/random_rotation/strided_slice_2/stack_1:output:0;sequential/random_rotation/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2,
*sequential/random_rotation/strided_slice_2│
!sequential/random_rotation/Cast_1Cast3sequential/random_rotation/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2#
!sequential/random_rotation/Cast_1╧
1sequential/random_rotation/stateful_uniform/shapePack1sequential/random_rotation/strided_slice:output:0*
N*
T0*
_output_shapes
:23
1sequential/random_rotation/stateful_uniform/shapeз
/sequential/random_rotation/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *█╔└21
/sequential/random_rotation/stateful_uniform/minз
/sequential/random_rotation/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *█╔@21
/sequential/random_rotation/stateful_uniform/max░
1sequential/random_rotation/stateful_uniform/ConstConst*
_output_shapes
:*
dtype0*
valueB: 23
1sequential/random_rotation/stateful_uniform/ConstЕ
0sequential/random_rotation/stateful_uniform/ProdProd:sequential/random_rotation/stateful_uniform/shape:output:0:sequential/random_rotation/stateful_uniform/Const:output:0*
T0*
_output_shapes
: 22
0sequential/random_rotation/stateful_uniform/Prodк
2sequential/random_rotation/stateful_uniform/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :24
2sequential/random_rotation/stateful_uniform/Cast/x█
2sequential/random_rotation/stateful_uniform/Cast_1Cast9sequential/random_rotation/stateful_uniform/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 24
2sequential/random_rotation/stateful_uniform/Cast_1р
:sequential/random_rotation/stateful_uniform/RngReadAndSkipRngReadAndSkipCsequential_random_rotation_stateful_uniform_rngreadandskip_resource;sequential/random_rotation/stateful_uniform/Cast/x:output:06sequential/random_rotation/stateful_uniform/Cast_1:y:0*
_output_shapes
:2<
:sequential/random_rotation/stateful_uniform/RngReadAndSkip╠
?sequential/random_rotation/stateful_uniform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2A
?sequential/random_rotation/stateful_uniform/strided_slice/stack╨
Asequential/random_rotation/stateful_uniform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2C
Asequential/random_rotation/stateful_uniform/strided_slice/stack_1╨
Asequential/random_rotation/stateful_uniform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2C
Asequential/random_rotation/stateful_uniform/strided_slice/stack_2Ё
9sequential/random_rotation/stateful_uniform/strided_sliceStridedSliceBsequential/random_rotation/stateful_uniform/RngReadAndSkip:value:0Hsequential/random_rotation/stateful_uniform/strided_slice/stack:output:0Jsequential/random_rotation/stateful_uniform/strided_slice/stack_1:output:0Jsequential/random_rotation/stateful_uniform/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2;
9sequential/random_rotation/stateful_uniform/strided_sliceъ
3sequential/random_rotation/stateful_uniform/BitcastBitcastBsequential/random_rotation/stateful_uniform/strided_slice:output:0*
T0	*
_output_shapes
:*

type025
3sequential/random_rotation/stateful_uniform/Bitcast╨
Asequential/random_rotation/stateful_uniform/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2C
Asequential/random_rotation/stateful_uniform/strided_slice_1/stack╘
Csequential/random_rotation/stateful_uniform/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2E
Csequential/random_rotation/stateful_uniform/strided_slice_1/stack_1╘
Csequential/random_rotation/stateful_uniform/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2E
Csequential/random_rotation/stateful_uniform/strided_slice_1/stack_2ш
;sequential/random_rotation/stateful_uniform/strided_slice_1StridedSliceBsequential/random_rotation/stateful_uniform/RngReadAndSkip:value:0Jsequential/random_rotation/stateful_uniform/strided_slice_1/stack:output:0Lsequential/random_rotation/stateful_uniform/strided_slice_1/stack_1:output:0Lsequential/random_rotation/stateful_uniform/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2=
;sequential/random_rotation/stateful_uniform/strided_slice_1Ё
5sequential/random_rotation/stateful_uniform/Bitcast_1BitcastDsequential/random_rotation/stateful_uniform/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type027
5sequential/random_rotation/stateful_uniform/Bitcast_1╓
Hsequential/random_rotation/stateful_uniform/StatelessRandomUniformV2/algConst*
_output_shapes
: *
dtype0*
value	B :2J
Hsequential/random_rotation/stateful_uniform/StatelessRandomUniformV2/alg┌
Dsequential/random_rotation/stateful_uniform/StatelessRandomUniformV2StatelessRandomUniformV2:sequential/random_rotation/stateful_uniform/shape:output:0>sequential/random_rotation/stateful_uniform/Bitcast_1:output:0<sequential/random_rotation/stateful_uniform/Bitcast:output:0Qsequential/random_rotation/stateful_uniform/StatelessRandomUniformV2/alg:output:0*#
_output_shapes
:         2F
Dsequential/random_rotation/stateful_uniform/StatelessRandomUniformV2■
/sequential/random_rotation/stateful_uniform/subSub8sequential/random_rotation/stateful_uniform/max:output:08sequential/random_rotation/stateful_uniform/min:output:0*
T0*
_output_shapes
: 21
/sequential/random_rotation/stateful_uniform/subЫ
/sequential/random_rotation/stateful_uniform/mulMulMsequential/random_rotation/stateful_uniform/StatelessRandomUniformV2:output:03sequential/random_rotation/stateful_uniform/sub:z:0*
T0*#
_output_shapes
:         21
/sequential/random_rotation/stateful_uniform/mulА
+sequential/random_rotation/stateful_uniformAddV23sequential/random_rotation/stateful_uniform/mul:z:08sequential/random_rotation/stateful_uniform/min:output:0*
T0*#
_output_shapes
:         2-
+sequential/random_rotation/stateful_uniformй
0sequential/random_rotation/rotation_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?22
0sequential/random_rotation/rotation_matrix/sub/yъ
.sequential/random_rotation/rotation_matrix/subSub%sequential/random_rotation/Cast_1:y:09sequential/random_rotation/rotation_matrix/sub/y:output:0*
T0*
_output_shapes
: 20
.sequential/random_rotation/rotation_matrix/sub╞
.sequential/random_rotation/rotation_matrix/CosCos/sequential/random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         20
.sequential/random_rotation/rotation_matrix/Cosн
2sequential/random_rotation/rotation_matrix/sub_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?24
2sequential/random_rotation/rotation_matrix/sub_1/yЁ
0sequential/random_rotation/rotation_matrix/sub_1Sub%sequential/random_rotation/Cast_1:y:0;sequential/random_rotation/rotation_matrix/sub_1/y:output:0*
T0*
_output_shapes
: 22
0sequential/random_rotation/rotation_matrix/sub_1 
.sequential/random_rotation/rotation_matrix/mulMul2sequential/random_rotation/rotation_matrix/Cos:y:04sequential/random_rotation/rotation_matrix/sub_1:z:0*
T0*#
_output_shapes
:         20
.sequential/random_rotation/rotation_matrix/mul╞
.sequential/random_rotation/rotation_matrix/SinSin/sequential/random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         20
.sequential/random_rotation/rotation_matrix/Sinн
2sequential/random_rotation/rotation_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?24
2sequential/random_rotation/rotation_matrix/sub_2/yю
0sequential/random_rotation/rotation_matrix/sub_2Sub#sequential/random_rotation/Cast:y:0;sequential/random_rotation/rotation_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 22
0sequential/random_rotation/rotation_matrix/sub_2Г
0sequential/random_rotation/rotation_matrix/mul_1Mul2sequential/random_rotation/rotation_matrix/Sin:y:04sequential/random_rotation/rotation_matrix/sub_2:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/mul_1Г
0sequential/random_rotation/rotation_matrix/sub_3Sub2sequential/random_rotation/rotation_matrix/mul:z:04sequential/random_rotation/rotation_matrix/mul_1:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/sub_3Г
0sequential/random_rotation/rotation_matrix/sub_4Sub2sequential/random_rotation/rotation_matrix/sub:z:04sequential/random_rotation/rotation_matrix/sub_3:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/sub_4▒
4sequential/random_rotation/rotation_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @26
4sequential/random_rotation/rotation_matrix/truediv/yЦ
2sequential/random_rotation/rotation_matrix/truedivRealDiv4sequential/random_rotation/rotation_matrix/sub_4:z:0=sequential/random_rotation/rotation_matrix/truediv/y:output:0*
T0*#
_output_shapes
:         24
2sequential/random_rotation/rotation_matrix/truedivн
2sequential/random_rotation/rotation_matrix/sub_5/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?24
2sequential/random_rotation/rotation_matrix/sub_5/yю
0sequential/random_rotation/rotation_matrix/sub_5Sub#sequential/random_rotation/Cast:y:0;sequential/random_rotation/rotation_matrix/sub_5/y:output:0*
T0*
_output_shapes
: 22
0sequential/random_rotation/rotation_matrix/sub_5╩
0sequential/random_rotation/rotation_matrix/Sin_1Sin/sequential/random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/Sin_1н
2sequential/random_rotation/rotation_matrix/sub_6/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?24
2sequential/random_rotation/rotation_matrix/sub_6/yЁ
0sequential/random_rotation/rotation_matrix/sub_6Sub%sequential/random_rotation/Cast_1:y:0;sequential/random_rotation/rotation_matrix/sub_6/y:output:0*
T0*
_output_shapes
: 22
0sequential/random_rotation/rotation_matrix/sub_6Е
0sequential/random_rotation/rotation_matrix/mul_2Mul4sequential/random_rotation/rotation_matrix/Sin_1:y:04sequential/random_rotation/rotation_matrix/sub_6:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/mul_2╩
0sequential/random_rotation/rotation_matrix/Cos_1Cos/sequential/random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/Cos_1н
2sequential/random_rotation/rotation_matrix/sub_7/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?24
2sequential/random_rotation/rotation_matrix/sub_7/yю
0sequential/random_rotation/rotation_matrix/sub_7Sub#sequential/random_rotation/Cast:y:0;sequential/random_rotation/rotation_matrix/sub_7/y:output:0*
T0*
_output_shapes
: 22
0sequential/random_rotation/rotation_matrix/sub_7Е
0sequential/random_rotation/rotation_matrix/mul_3Mul4sequential/random_rotation/rotation_matrix/Cos_1:y:04sequential/random_rotation/rotation_matrix/sub_7:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/mul_3Г
.sequential/random_rotation/rotation_matrix/addAddV24sequential/random_rotation/rotation_matrix/mul_2:z:04sequential/random_rotation/rotation_matrix/mul_3:z:0*
T0*#
_output_shapes
:         20
.sequential/random_rotation/rotation_matrix/addГ
0sequential/random_rotation/rotation_matrix/sub_8Sub4sequential/random_rotation/rotation_matrix/sub_5:z:02sequential/random_rotation/rotation_matrix/add:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/sub_8╡
6sequential/random_rotation/rotation_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @28
6sequential/random_rotation/rotation_matrix/truediv_1/yЬ
4sequential/random_rotation/rotation_matrix/truediv_1RealDiv4sequential/random_rotation/rotation_matrix/sub_8:z:0?sequential/random_rotation/rotation_matrix/truediv_1/y:output:0*
T0*#
_output_shapes
:         26
4sequential/random_rotation/rotation_matrix/truediv_1├
0sequential/random_rotation/rotation_matrix/ShapeShape/sequential/random_rotation/stateful_uniform:z:0*
T0*
_output_shapes
:22
0sequential/random_rotation/rotation_matrix/Shape╩
>sequential/random_rotation/rotation_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2@
>sequential/random_rotation/rotation_matrix/strided_slice/stack╬
@sequential/random_rotation/rotation_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2B
@sequential/random_rotation/rotation_matrix/strided_slice/stack_1╬
@sequential/random_rotation/rotation_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2B
@sequential/random_rotation/rotation_matrix/strided_slice/stack_2ф
8sequential/random_rotation/rotation_matrix/strided_sliceStridedSlice9sequential/random_rotation/rotation_matrix/Shape:output:0Gsequential/random_rotation/rotation_matrix/strided_slice/stack:output:0Isequential/random_rotation/rotation_matrix/strided_slice/stack_1:output:0Isequential/random_rotation/rotation_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2:
8sequential/random_rotation/rotation_matrix/strided_slice╩
0sequential/random_rotation/rotation_matrix/Cos_2Cos/sequential/random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/Cos_2╒
@sequential/random_rotation/rotation_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB"        2B
@sequential/random_rotation/rotation_matrix/strided_slice_1/stack┘
Bsequential/random_rotation/rotation_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2D
Bsequential/random_rotation/rotation_matrix/strided_slice_1/stack_1┘
Bsequential/random_rotation/rotation_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2D
Bsequential/random_rotation/rotation_matrix/strided_slice_1/stack_2Щ
:sequential/random_rotation/rotation_matrix/strided_slice_1StridedSlice4sequential/random_rotation/rotation_matrix/Cos_2:y:0Isequential/random_rotation/rotation_matrix/strided_slice_1/stack:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_1/stack_1:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2<
:sequential/random_rotation/rotation_matrix/strided_slice_1╩
0sequential/random_rotation/rotation_matrix/Sin_2Sin/sequential/random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/Sin_2╒
@sequential/random_rotation/rotation_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB"        2B
@sequential/random_rotation/rotation_matrix/strided_slice_2/stack┘
Bsequential/random_rotation/rotation_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2D
Bsequential/random_rotation/rotation_matrix/strided_slice_2/stack_1┘
Bsequential/random_rotation/rotation_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2D
Bsequential/random_rotation/rotation_matrix/strided_slice_2/stack_2Щ
:sequential/random_rotation/rotation_matrix/strided_slice_2StridedSlice4sequential/random_rotation/rotation_matrix/Sin_2:y:0Isequential/random_rotation/rotation_matrix/strided_slice_2/stack:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_2/stack_1:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2<
:sequential/random_rotation/rotation_matrix/strided_slice_2▐
.sequential/random_rotation/rotation_matrix/NegNegCsequential/random_rotation/rotation_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:         20
.sequential/random_rotation/rotation_matrix/Neg╒
@sequential/random_rotation/rotation_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB"        2B
@sequential/random_rotation/rotation_matrix/strided_slice_3/stack┘
Bsequential/random_rotation/rotation_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2D
Bsequential/random_rotation/rotation_matrix/strided_slice_3/stack_1┘
Bsequential/random_rotation/rotation_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2D
Bsequential/random_rotation/rotation_matrix/strided_slice_3/stack_2Ы
:sequential/random_rotation/rotation_matrix/strided_slice_3StridedSlice6sequential/random_rotation/rotation_matrix/truediv:z:0Isequential/random_rotation/rotation_matrix/strided_slice_3/stack:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_3/stack_1:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2<
:sequential/random_rotation/rotation_matrix/strided_slice_3╩
0sequential/random_rotation/rotation_matrix/Sin_3Sin/sequential/random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/Sin_3╒
@sequential/random_rotation/rotation_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*
valueB"        2B
@sequential/random_rotation/rotation_matrix/strided_slice_4/stack┘
Bsequential/random_rotation/rotation_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2D
Bsequential/random_rotation/rotation_matrix/strided_slice_4/stack_1┘
Bsequential/random_rotation/rotation_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2D
Bsequential/random_rotation/rotation_matrix/strided_slice_4/stack_2Щ
:sequential/random_rotation/rotation_matrix/strided_slice_4StridedSlice4sequential/random_rotation/rotation_matrix/Sin_3:y:0Isequential/random_rotation/rotation_matrix/strided_slice_4/stack:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_4/stack_1:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2<
:sequential/random_rotation/rotation_matrix/strided_slice_4╩
0sequential/random_rotation/rotation_matrix/Cos_3Cos/sequential/random_rotation/stateful_uniform:z:0*
T0*#
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/Cos_3╒
@sequential/random_rotation/rotation_matrix/strided_slice_5/stackConst*
_output_shapes
:*
dtype0*
valueB"        2B
@sequential/random_rotation/rotation_matrix/strided_slice_5/stack┘
Bsequential/random_rotation/rotation_matrix/strided_slice_5/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2D
Bsequential/random_rotation/rotation_matrix/strided_slice_5/stack_1┘
Bsequential/random_rotation/rotation_matrix/strided_slice_5/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2D
Bsequential/random_rotation/rotation_matrix/strided_slice_5/stack_2Щ
:sequential/random_rotation/rotation_matrix/strided_slice_5StridedSlice4sequential/random_rotation/rotation_matrix/Cos_3:y:0Isequential/random_rotation/rotation_matrix/strided_slice_5/stack:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_5/stack_1:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_5/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2<
:sequential/random_rotation/rotation_matrix/strided_slice_5╒
@sequential/random_rotation/rotation_matrix/strided_slice_6/stackConst*
_output_shapes
:*
dtype0*
valueB"        2B
@sequential/random_rotation/rotation_matrix/strided_slice_6/stack┘
Bsequential/random_rotation/rotation_matrix/strided_slice_6/stack_1Const*
_output_shapes
:*
dtype0*
valueB"        2D
Bsequential/random_rotation/rotation_matrix/strided_slice_6/stack_1┘
Bsequential/random_rotation/rotation_matrix/strided_slice_6/stack_2Const*
_output_shapes
:*
dtype0*
valueB"      2D
Bsequential/random_rotation/rotation_matrix/strided_slice_6/stack_2Э
:sequential/random_rotation/rotation_matrix/strided_slice_6StridedSlice8sequential/random_rotation/rotation_matrix/truediv_1:z:0Isequential/random_rotation/rotation_matrix/strided_slice_6/stack:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_6/stack_1:output:0Ksequential/random_rotation/rotation_matrix/strided_slice_6/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask2<
:sequential/random_rotation/rotation_matrix/strided_slice_6▓
6sequential/random_rotation/rotation_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :28
6sequential/random_rotation/rotation_matrix/zeros/mul/yШ
4sequential/random_rotation/rotation_matrix/zeros/mulMulAsequential/random_rotation/rotation_matrix/strided_slice:output:0?sequential/random_rotation/rotation_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 26
4sequential/random_rotation/rotation_matrix/zeros/mul╡
7sequential/random_rotation/rotation_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш29
7sequential/random_rotation/rotation_matrix/zeros/Less/yУ
5sequential/random_rotation/rotation_matrix/zeros/LessLess8sequential/random_rotation/rotation_matrix/zeros/mul:z:0@sequential/random_rotation/rotation_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 27
5sequential/random_rotation/rotation_matrix/zeros/Less╕
9sequential/random_rotation/rotation_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2;
9sequential/random_rotation/rotation_matrix/zeros/packed/1п
7sequential/random_rotation/rotation_matrix/zeros/packedPackAsequential/random_rotation/rotation_matrix/strided_slice:output:0Bsequential/random_rotation/rotation_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:29
7sequential/random_rotation/rotation_matrix/zeros/packed╡
6sequential/random_rotation/rotation_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    28
6sequential/random_rotation/rotation_matrix/zeros/Constб
0sequential/random_rotation/rotation_matrix/zerosFill@sequential/random_rotation/rotation_matrix/zeros/packed:output:0?sequential/random_rotation/rotation_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:         22
0sequential/random_rotation/rotation_matrix/zeros▓
6sequential/random_rotation/rotation_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :28
6sequential/random_rotation/rotation_matrix/concat/axis╢
1sequential/random_rotation/rotation_matrix/concatConcatV2Csequential/random_rotation/rotation_matrix/strided_slice_1:output:02sequential/random_rotation/rotation_matrix/Neg:y:0Csequential/random_rotation/rotation_matrix/strided_slice_3:output:0Csequential/random_rotation/rotation_matrix/strided_slice_4:output:0Csequential/random_rotation/rotation_matrix/strided_slice_5:output:0Csequential/random_rotation/rotation_matrix/strided_slice_6:output:09sequential/random_rotation/rotation_matrix/zeros:output:0?sequential/random_rotation/rotation_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:         23
1sequential/random_rotation/rotation_matrix/concat╟
*sequential/random_rotation/transform/ShapeShape?sequential/random_flip/stateless_random_flip_left_right/add:z:0*
T0*
_output_shapes
:2,
*sequential/random_rotation/transform/Shape╛
8sequential/random_rotation/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2:
8sequential/random_rotation/transform/strided_slice/stack┬
:sequential/random_rotation/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2<
:sequential/random_rotation/transform/strided_slice/stack_1┬
:sequential/random_rotation/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2<
:sequential/random_rotation/transform/strided_slice/stack_2м
2sequential/random_rotation/transform/strided_sliceStridedSlice3sequential/random_rotation/transform/Shape:output:0Asequential/random_rotation/transform/strided_slice/stack:output:0Csequential/random_rotation/transform/strided_slice/stack_1:output:0Csequential/random_rotation/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:24
2sequential/random_rotation/transform/strided_sliceз
/sequential/random_rotation/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    21
/sequential/random_rotation/transform/fill_valueЗ
?sequential/random_rotation/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3?sequential/random_flip/stateless_random_flip_left_right/add:z:0:sequential/random_rotation/rotation_matrix/concat:output:0;sequential/random_rotation/transform/strided_slice:output:08sequential/random_rotation/transform/fill_value:output:0*/
_output_shapes
:         &*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2A
?sequential/random_rotation/transform/ImageProjectiveTransformV3└
sequential/random_zoom/ShapeShapeTsequential/random_rotation/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2
sequential/random_zoom/Shapeв
*sequential/random_zoom/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2,
*sequential/random_zoom/strided_slice/stackж
,sequential/random_zoom/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2.
,sequential/random_zoom/strided_slice/stack_1ж
,sequential/random_zoom/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,sequential/random_zoom/strided_slice/stack_2ь
$sequential/random_zoom/strided_sliceStridedSlice%sequential/random_zoom/Shape:output:03sequential/random_zoom/strided_slice/stack:output:05sequential/random_zoom/strided_slice/stack_1:output:05sequential/random_zoom/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2&
$sequential/random_zoom/strided_sliceп
,sequential/random_zoom/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:
¤        2.
,sequential/random_zoom/strided_slice_1/stack│
.sequential/random_zoom/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
■        20
.sequential/random_zoom/strided_slice_1/stack_1к
.sequential/random_zoom/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:20
.sequential/random_zoom/strided_slice_1/stack_2Ў
&sequential/random_zoom/strided_slice_1StridedSlice%sequential/random_zoom/Shape:output:05sequential/random_zoom/strided_slice_1/stack:output:07sequential/random_zoom/strided_slice_1/stack_1:output:07sequential/random_zoom/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2(
&sequential/random_zoom/strided_slice_1г
sequential/random_zoom/CastCast/sequential/random_zoom/strided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
sequential/random_zoom/Castп
,sequential/random_zoom/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2.
,sequential/random_zoom/strided_slice_2/stack│
.sequential/random_zoom/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         20
.sequential/random_zoom/strided_slice_2/stack_1к
.sequential/random_zoom/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:20
.sequential/random_zoom/strided_slice_2/stack_2Ў
&sequential/random_zoom/strided_slice_2StridedSlice%sequential/random_zoom/Shape:output:05sequential/random_zoom/strided_slice_2/stack:output:07sequential/random_zoom/strided_slice_2/stack_1:output:07sequential/random_zoom/strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2(
&sequential/random_zoom/strided_slice_2з
sequential/random_zoom/Cast_1Cast/sequential/random_zoom/strided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
sequential/random_zoom/Cast_1д
/sequential/random_zoom/stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :21
/sequential/random_zoom/stateful_uniform/shape/1¤
-sequential/random_zoom/stateful_uniform/shapePack-sequential/random_zoom/strided_slice:output:08sequential/random_zoom/stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2/
-sequential/random_zoom/stateful_uniform/shapeЯ
+sequential/random_zoom/stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *fff?2-
+sequential/random_zoom/stateful_uniform/minЯ
+sequential/random_zoom/stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *═╠М?2-
+sequential/random_zoom/stateful_uniform/maxи
-sequential/random_zoom/stateful_uniform/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2/
-sequential/random_zoom/stateful_uniform/Constї
,sequential/random_zoom/stateful_uniform/ProdProd6sequential/random_zoom/stateful_uniform/shape:output:06sequential/random_zoom/stateful_uniform/Const:output:0*
T0*
_output_shapes
: 2.
,sequential/random_zoom/stateful_uniform/Prodв
.sequential/random_zoom/stateful_uniform/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :20
.sequential/random_zoom/stateful_uniform/Cast/x╧
.sequential/random_zoom/stateful_uniform/Cast_1Cast5sequential/random_zoom/stateful_uniform/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 20
.sequential/random_zoom/stateful_uniform/Cast_1╠
6sequential/random_zoom/stateful_uniform/RngReadAndSkipRngReadAndSkip?sequential_random_zoom_stateful_uniform_rngreadandskip_resource7sequential/random_zoom/stateful_uniform/Cast/x:output:02sequential/random_zoom/stateful_uniform/Cast_1:y:0*
_output_shapes
:28
6sequential/random_zoom/stateful_uniform/RngReadAndSkip─
;sequential/random_zoom/stateful_uniform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2=
;sequential/random_zoom/stateful_uniform/strided_slice/stack╚
=sequential/random_zoom/stateful_uniform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2?
=sequential/random_zoom/stateful_uniform/strided_slice/stack_1╚
=sequential/random_zoom/stateful_uniform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2?
=sequential/random_zoom/stateful_uniform/strided_slice/stack_2╪
5sequential/random_zoom/stateful_uniform/strided_sliceStridedSlice>sequential/random_zoom/stateful_uniform/RngReadAndSkip:value:0Dsequential/random_zoom/stateful_uniform/strided_slice/stack:output:0Fsequential/random_zoom/stateful_uniform/strided_slice/stack_1:output:0Fsequential/random_zoom/stateful_uniform/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask27
5sequential/random_zoom/stateful_uniform/strided_slice▐
/sequential/random_zoom/stateful_uniform/BitcastBitcast>sequential/random_zoom/stateful_uniform/strided_slice:output:0*
T0	*
_output_shapes
:*

type021
/sequential/random_zoom/stateful_uniform/Bitcast╚
=sequential/random_zoom/stateful_uniform/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2?
=sequential/random_zoom/stateful_uniform/strided_slice_1/stack╠
?sequential/random_zoom/stateful_uniform/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2A
?sequential/random_zoom/stateful_uniform/strided_slice_1/stack_1╠
?sequential/random_zoom/stateful_uniform/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2A
?sequential/random_zoom/stateful_uniform/strided_slice_1/stack_2╨
7sequential/random_zoom/stateful_uniform/strided_slice_1StridedSlice>sequential/random_zoom/stateful_uniform/RngReadAndSkip:value:0Fsequential/random_zoom/stateful_uniform/strided_slice_1/stack:output:0Hsequential/random_zoom/stateful_uniform/strided_slice_1/stack_1:output:0Hsequential/random_zoom/stateful_uniform/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:29
7sequential/random_zoom/stateful_uniform/strided_slice_1ф
1sequential/random_zoom/stateful_uniform/Bitcast_1Bitcast@sequential/random_zoom/stateful_uniform/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type023
1sequential/random_zoom/stateful_uniform/Bitcast_1╬
Dsequential/random_zoom/stateful_uniform/StatelessRandomUniformV2/algConst*
_output_shapes
: *
dtype0*
value	B :2F
Dsequential/random_zoom/stateful_uniform/StatelessRandomUniformV2/alg╞
@sequential/random_zoom/stateful_uniform/StatelessRandomUniformV2StatelessRandomUniformV26sequential/random_zoom/stateful_uniform/shape:output:0:sequential/random_zoom/stateful_uniform/Bitcast_1:output:08sequential/random_zoom/stateful_uniform/Bitcast:output:0Msequential/random_zoom/stateful_uniform/StatelessRandomUniformV2/alg:output:0*'
_output_shapes
:         2B
@sequential/random_zoom/stateful_uniform/StatelessRandomUniformV2ю
+sequential/random_zoom/stateful_uniform/subSub4sequential/random_zoom/stateful_uniform/max:output:04sequential/random_zoom/stateful_uniform/min:output:0*
T0*
_output_shapes
: 2-
+sequential/random_zoom/stateful_uniform/subП
+sequential/random_zoom/stateful_uniform/mulMulIsequential/random_zoom/stateful_uniform/StatelessRandomUniformV2:output:0/sequential/random_zoom/stateful_uniform/sub:z:0*
T0*'
_output_shapes
:         2-
+sequential/random_zoom/stateful_uniform/mulЇ
'sequential/random_zoom/stateful_uniformAddV2/sequential/random_zoom/stateful_uniform/mul:z:04sequential/random_zoom/stateful_uniform/min:output:0*
T0*'
_output_shapes
:         2)
'sequential/random_zoom/stateful_uniformК
"sequential/random_zoom/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2$
"sequential/random_zoom/concat/axisМ
sequential/random_zoom/concatConcatV2+sequential/random_zoom/stateful_uniform:z:0+sequential/random_zoom/stateful_uniform:z:0+sequential/random_zoom/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2
sequential/random_zoom/concatк
(sequential/random_zoom/zoom_matrix/ShapeShape&sequential/random_zoom/concat:output:0*
T0*
_output_shapes
:2*
(sequential/random_zoom/zoom_matrix/Shape║
6sequential/random_zoom/zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 28
6sequential/random_zoom/zoom_matrix/strided_slice/stack╛
8sequential/random_zoom/zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2:
8sequential/random_zoom/zoom_matrix/strided_slice/stack_1╛
8sequential/random_zoom/zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2:
8sequential/random_zoom/zoom_matrix/strided_slice/stack_2┤
0sequential/random_zoom/zoom_matrix/strided_sliceStridedSlice1sequential/random_zoom/zoom_matrix/Shape:output:0?sequential/random_zoom/zoom_matrix/strided_slice/stack:output:0Asequential/random_zoom/zoom_matrix/strided_slice/stack_1:output:0Asequential/random_zoom/zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask22
0sequential/random_zoom/zoom_matrix/strided_sliceЩ
(sequential/random_zoom/zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2*
(sequential/random_zoom/zoom_matrix/sub/y╬
&sequential/random_zoom/zoom_matrix/subSub!sequential/random_zoom/Cast_1:y:01sequential/random_zoom/zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 2(
&sequential/random_zoom/zoom_matrix/subб
,sequential/random_zoom/zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2.
,sequential/random_zoom/zoom_matrix/truediv/yч
*sequential/random_zoom/zoom_matrix/truedivRealDiv*sequential/random_zoom/zoom_matrix/sub:z:05sequential/random_zoom/zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 2,
*sequential/random_zoom/zoom_matrix/truediv╔
8sequential/random_zoom/zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2:
8sequential/random_zoom/zoom_matrix/strided_slice_1/stack═
:sequential/random_zoom/zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2<
:sequential/random_zoom/zoom_matrix/strided_slice_1/stack_1═
:sequential/random_zoom/zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2<
:sequential/random_zoom/zoom_matrix/strided_slice_1/stack_2√
2sequential/random_zoom/zoom_matrix/strided_slice_1StridedSlice&sequential/random_zoom/concat:output:0Asequential/random_zoom/zoom_matrix/strided_slice_1/stack:output:0Csequential/random_zoom/zoom_matrix/strided_slice_1/stack_1:output:0Csequential/random_zoom/zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask24
2sequential/random_zoom/zoom_matrix/strided_slice_1Э
*sequential/random_zoom/zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2,
*sequential/random_zoom/zoom_matrix/sub_1/x 
(sequential/random_zoom/zoom_matrix/sub_1Sub3sequential/random_zoom/zoom_matrix/sub_1/x:output:0;sequential/random_zoom/zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:         2*
(sequential/random_zoom/zoom_matrix/sub_1ч
&sequential/random_zoom/zoom_matrix/mulMul.sequential/random_zoom/zoom_matrix/truediv:z:0,sequential/random_zoom/zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:         2(
&sequential/random_zoom/zoom_matrix/mulЭ
*sequential/random_zoom/zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2,
*sequential/random_zoom/zoom_matrix/sub_2/y╥
(sequential/random_zoom/zoom_matrix/sub_2Subsequential/random_zoom/Cast:y:03sequential/random_zoom/zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2*
(sequential/random_zoom/zoom_matrix/sub_2е
.sequential/random_zoom/zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @20
.sequential/random_zoom/zoom_matrix/truediv_1/yя
,sequential/random_zoom/zoom_matrix/truediv_1RealDiv,sequential/random_zoom/zoom_matrix/sub_2:z:07sequential/random_zoom/zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 2.
,sequential/random_zoom/zoom_matrix/truediv_1╔
8sequential/random_zoom/zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2:
8sequential/random_zoom/zoom_matrix/strided_slice_2/stack═
:sequential/random_zoom/zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2<
:sequential/random_zoom/zoom_matrix/strided_slice_2/stack_1═
:sequential/random_zoom/zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2<
:sequential/random_zoom/zoom_matrix/strided_slice_2/stack_2√
2sequential/random_zoom/zoom_matrix/strided_slice_2StridedSlice&sequential/random_zoom/concat:output:0Asequential/random_zoom/zoom_matrix/strided_slice_2/stack:output:0Csequential/random_zoom/zoom_matrix/strided_slice_2/stack_1:output:0Csequential/random_zoom/zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask24
2sequential/random_zoom/zoom_matrix/strided_slice_2Э
*sequential/random_zoom/zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2,
*sequential/random_zoom/zoom_matrix/sub_3/x 
(sequential/random_zoom/zoom_matrix/sub_3Sub3sequential/random_zoom/zoom_matrix/sub_3/x:output:0;sequential/random_zoom/zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:         2*
(sequential/random_zoom/zoom_matrix/sub_3э
(sequential/random_zoom/zoom_matrix/mul_1Mul0sequential/random_zoom/zoom_matrix/truediv_1:z:0,sequential/random_zoom/zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:         2*
(sequential/random_zoom/zoom_matrix/mul_1╔
8sequential/random_zoom/zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2:
8sequential/random_zoom/zoom_matrix/strided_slice_3/stack═
:sequential/random_zoom/zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2<
:sequential/random_zoom/zoom_matrix/strided_slice_3/stack_1═
:sequential/random_zoom/zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2<
:sequential/random_zoom/zoom_matrix/strided_slice_3/stack_2√
2sequential/random_zoom/zoom_matrix/strided_slice_3StridedSlice&sequential/random_zoom/concat:output:0Asequential/random_zoom/zoom_matrix/strided_slice_3/stack:output:0Csequential/random_zoom/zoom_matrix/strided_slice_3/stack_1:output:0Csequential/random_zoom/zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask24
2sequential/random_zoom/zoom_matrix/strided_slice_3в
.sequential/random_zoom/zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :20
.sequential/random_zoom/zoom_matrix/zeros/mul/y°
,sequential/random_zoom/zoom_matrix/zeros/mulMul9sequential/random_zoom/zoom_matrix/strided_slice:output:07sequential/random_zoom/zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2.
,sequential/random_zoom/zoom_matrix/zeros/mulе
/sequential/random_zoom/zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш21
/sequential/random_zoom/zoom_matrix/zeros/Less/yє
-sequential/random_zoom/zoom_matrix/zeros/LessLess0sequential/random_zoom/zoom_matrix/zeros/mul:z:08sequential/random_zoom/zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2/
-sequential/random_zoom/zoom_matrix/zeros/Lessи
1sequential/random_zoom/zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :23
1sequential/random_zoom/zoom_matrix/zeros/packed/1П
/sequential/random_zoom/zoom_matrix/zeros/packedPack9sequential/random_zoom/zoom_matrix/strided_slice:output:0:sequential/random_zoom/zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:21
/sequential/random_zoom/zoom_matrix/zeros/packedе
.sequential/random_zoom/zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    20
.sequential/random_zoom/zoom_matrix/zeros/ConstБ
(sequential/random_zoom/zoom_matrix/zerosFill8sequential/random_zoom/zoom_matrix/zeros/packed:output:07sequential/random_zoom/zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:         2*
(sequential/random_zoom/zoom_matrix/zerosж
0sequential/random_zoom/zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :22
0sequential/random_zoom/zoom_matrix/zeros_1/mul/y■
.sequential/random_zoom/zoom_matrix/zeros_1/mulMul9sequential/random_zoom/zoom_matrix/strided_slice:output:09sequential/random_zoom/zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 20
.sequential/random_zoom/zoom_matrix/zeros_1/mulй
1sequential/random_zoom/zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш23
1sequential/random_zoom/zoom_matrix/zeros_1/Less/y√
/sequential/random_zoom/zoom_matrix/zeros_1/LessLess2sequential/random_zoom/zoom_matrix/zeros_1/mul:z:0:sequential/random_zoom/zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 21
/sequential/random_zoom/zoom_matrix/zeros_1/Lessм
3sequential/random_zoom/zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :25
3sequential/random_zoom/zoom_matrix/zeros_1/packed/1Х
1sequential/random_zoom/zoom_matrix/zeros_1/packedPack9sequential/random_zoom/zoom_matrix/strided_slice:output:0<sequential/random_zoom/zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:23
1sequential/random_zoom/zoom_matrix/zeros_1/packedй
0sequential/random_zoom/zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    22
0sequential/random_zoom/zoom_matrix/zeros_1/ConstЙ
*sequential/random_zoom/zoom_matrix/zeros_1Fill:sequential/random_zoom/zoom_matrix/zeros_1/packed:output:09sequential/random_zoom/zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:         2,
*sequential/random_zoom/zoom_matrix/zeros_1╔
8sequential/random_zoom/zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2:
8sequential/random_zoom/zoom_matrix/strided_slice_4/stack═
:sequential/random_zoom/zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2<
:sequential/random_zoom/zoom_matrix/strided_slice_4/stack_1═
:sequential/random_zoom/zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2<
:sequential/random_zoom/zoom_matrix/strided_slice_4/stack_2√
2sequential/random_zoom/zoom_matrix/strided_slice_4StridedSlice&sequential/random_zoom/concat:output:0Asequential/random_zoom/zoom_matrix/strided_slice_4/stack:output:0Csequential/random_zoom/zoom_matrix/strided_slice_4/stack_1:output:0Csequential/random_zoom/zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask24
2sequential/random_zoom/zoom_matrix/strided_slice_4ж
0sequential/random_zoom/zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :22
0sequential/random_zoom/zoom_matrix/zeros_2/mul/y■
.sequential/random_zoom/zoom_matrix/zeros_2/mulMul9sequential/random_zoom/zoom_matrix/strided_slice:output:09sequential/random_zoom/zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 20
.sequential/random_zoom/zoom_matrix/zeros_2/mulй
1sequential/random_zoom/zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш23
1sequential/random_zoom/zoom_matrix/zeros_2/Less/y√
/sequential/random_zoom/zoom_matrix/zeros_2/LessLess2sequential/random_zoom/zoom_matrix/zeros_2/mul:z:0:sequential/random_zoom/zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 21
/sequential/random_zoom/zoom_matrix/zeros_2/Lessм
3sequential/random_zoom/zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :25
3sequential/random_zoom/zoom_matrix/zeros_2/packed/1Х
1sequential/random_zoom/zoom_matrix/zeros_2/packedPack9sequential/random_zoom/zoom_matrix/strided_slice:output:0<sequential/random_zoom/zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:23
1sequential/random_zoom/zoom_matrix/zeros_2/packedй
0sequential/random_zoom/zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    22
0sequential/random_zoom/zoom_matrix/zeros_2/ConstЙ
*sequential/random_zoom/zoom_matrix/zeros_2Fill:sequential/random_zoom/zoom_matrix/zeros_2/packed:output:09sequential/random_zoom/zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:         2,
*sequential/random_zoom/zoom_matrix/zeros_2в
.sequential/random_zoom/zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :20
.sequential/random_zoom/zoom_matrix/concat/axis╟
)sequential/random_zoom/zoom_matrix/concatConcatV2;sequential/random_zoom/zoom_matrix/strided_slice_3:output:01sequential/random_zoom/zoom_matrix/zeros:output:0*sequential/random_zoom/zoom_matrix/mul:z:03sequential/random_zoom/zoom_matrix/zeros_1:output:0;sequential/random_zoom/zoom_matrix/strided_slice_4:output:0,sequential/random_zoom/zoom_matrix/mul_1:z:03sequential/random_zoom/zoom_matrix/zeros_2:output:07sequential/random_zoom/zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2+
)sequential/random_zoom/zoom_matrix/concat╘
&sequential/random_zoom/transform/ShapeShapeTsequential/random_rotation/transform/ImageProjectiveTransformV3:transformed_images:0*
T0*
_output_shapes
:2(
&sequential/random_zoom/transform/Shape╢
4sequential/random_zoom/transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:26
4sequential/random_zoom/transform/strided_slice/stack║
6sequential/random_zoom/transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:28
6sequential/random_zoom/transform/strided_slice/stack_1║
6sequential/random_zoom/transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:28
6sequential/random_zoom/transform/strided_slice/stack_2Ф
.sequential/random_zoom/transform/strided_sliceStridedSlice/sequential/random_zoom/transform/Shape:output:0=sequential/random_zoom/transform/strided_slice/stack:output:0?sequential/random_zoom/transform/strided_slice/stack_1:output:0?sequential/random_zoom/transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:20
.sequential/random_zoom/transform/strided_sliceЯ
+sequential/random_zoom/transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2-
+sequential/random_zoom/transform/fill_valueД
;sequential/random_zoom/transform/ImageProjectiveTransformV3ImageProjectiveTransformV3Tsequential/random_rotation/transform/ImageProjectiveTransformV3:transformed_images:02sequential/random_zoom/zoom_matrix/concat:output:07sequential/random_zoom/transform/strided_slice:output:04sequential/random_zoom/transform/fill_value:output:0*/
_output_shapes
:         &*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2=
;sequential/random_zoom/transform/ImageProjectiveTransformV3m
rescaling_1/Cast/xConst*
_output_shapes
: *
dtype0*
valueB
 *БАА;2
rescaling_1/Cast/xq
rescaling_1/Cast_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2
rescaling_1/Cast_1/x╥
rescaling_1/mulMulPsequential/random_zoom/transform/ImageProjectiveTransformV3:transformed_images:0rescaling_1/Cast/x:output:0*
T0*/
_output_shapes
:         &2
rescaling_1/mulЩ
rescaling_1/addAddV2rescaling_1/mul:z:0rescaling_1/Cast_1/x:output:0*
T0*/
_output_shapes
:         &2
rescaling_1/addк
conv2d/Conv2D/ReadVariableOpReadVariableOp%conv2d_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
conv2d/Conv2D/ReadVariableOp┼
conv2d/Conv2DConv2Drescaling_1/add:z:0$conv2d/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &*
paddingSAME*
strides
2
conv2d/Conv2Dб
conv2d/BiasAdd/ReadVariableOpReadVariableOp&conv2d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
conv2d/BiasAdd/ReadVariableOpд
conv2d/BiasAddBiasAddconv2d/Conv2D:output:0%conv2d/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &2
conv2d/BiasAddu
conv2d/ReluReluconv2d/BiasAdd:output:0*
T0*/
_output_shapes
:         &2
conv2d/Relu┴
max_pooling2d/MaxPoolMaxPoolconv2d/Relu:activations:0*/
_output_shapes
:         *
ksize
*
paddingVALID*
strides
2
max_pooling2d/MaxPool░
conv2d_1/Conv2D/ReadVariableOpReadVariableOp'conv2d_1_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02 
conv2d_1/Conv2D/ReadVariableOp╓
conv2d_1/Conv2DConv2Dmax_pooling2d/MaxPool:output:0&conv2d_1/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         *
paddingSAME*
strides
2
conv2d_1/Conv2Dз
conv2d_1/BiasAdd/ReadVariableOpReadVariableOp(conv2d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv2d_1/BiasAdd/ReadVariableOpм
conv2d_1/BiasAddBiasAddconv2d_1/Conv2D:output:0'conv2d_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
conv2d_1/BiasAdd{
conv2d_1/ReluReluconv2d_1/BiasAdd:output:0*
T0*/
_output_shapes
:         2
conv2d_1/Relu╟
max_pooling2d_1/MaxPoolMaxPoolconv2d_1/Relu:activations:0*/
_output_shapes
:         	*
ksize
*
paddingVALID*
strides
2
max_pooling2d_1/MaxPool░
conv2d_2/Conv2D/ReadVariableOpReadVariableOp'conv2d_2_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype02 
conv2d_2/Conv2D/ReadVariableOp╪
conv2d_2/Conv2DConv2D max_pooling2d_1/MaxPool:output:0&conv2d_2/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 *
paddingSAME*
strides
2
conv2d_2/Conv2Dз
conv2d_2/BiasAdd/ReadVariableOpReadVariableOp(conv2d_2_biasadd_readvariableop_resource*
_output_shapes
: *
dtype02!
conv2d_2/BiasAdd/ReadVariableOpм
conv2d_2/BiasAddBiasAddconv2d_2/Conv2D:output:0'conv2d_2/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 2
conv2d_2/BiasAdd{
conv2d_2/ReluReluconv2d_2/BiasAdd:output:0*
T0*/
_output_shapes
:         	 2
conv2d_2/Relu╟
max_pooling2d_2/MaxPoolMaxPoolconv2d_2/Relu:activations:0*/
_output_shapes
:          *
ksize
*
paddingVALID*
strides
2
max_pooling2d_2/MaxPool░
conv2d_3/Conv2D/ReadVariableOpReadVariableOp'conv2d_3_conv2d_readvariableop_resource*&
_output_shapes
: @*
dtype02 
conv2d_3/Conv2D/ReadVariableOp╪
conv2d_3/Conv2DConv2D max_pooling2d_2/MaxPool:output:0&conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @*
paddingSAME*
strides
2
conv2d_3/Conv2Dз
conv2d_3/BiasAdd/ReadVariableOpReadVariableOp(conv2d_3_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02!
conv2d_3/BiasAdd/ReadVariableOpм
conv2d_3/BiasAddBiasAddconv2d_3/Conv2D:output:0'conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         @2
conv2d_3/BiasAdd{
conv2d_3/ReluReluconv2d_3/BiasAdd:output:0*
T0*/
_output_shapes
:         @2
conv2d_3/Relu░
conv2d_4/Conv2D/ReadVariableOpReadVariableOp'conv2d_4_conv2d_readvariableop_resource*&
_output_shapes
:@P*
dtype02 
conv2d_4/Conv2D/ReadVariableOp╙
conv2d_4/Conv2DConv2Dconv2d_3/Relu:activations:0&conv2d_4/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P*
paddingSAME*
strides
2
conv2d_4/Conv2Dз
conv2d_4/BiasAdd/ReadVariableOpReadVariableOp(conv2d_4_biasadd_readvariableop_resource*
_output_shapes
:P*
dtype02!
conv2d_4/BiasAdd/ReadVariableOpм
conv2d_4/BiasAddBiasAddconv2d_4/Conv2D:output:0'conv2d_4/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         P2
conv2d_4/BiasAdd{
conv2d_4/ReluReluconv2d_4/BiasAdd:output:0*
T0*/
_output_shapes
:         P2
conv2d_4/Relu╟
max_pooling2d_3/MaxPoolMaxPoolconv2d_4/Relu:activations:0*/
_output_shapes
:         P*
ksize
*
paddingVALID*
strides
2
max_pooling2d_3/MaxPoolo
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"    а   2
flatten/ConstЪ
flatten/ReshapeReshape max_pooling2d_3/MaxPool:output:0flatten/Const:output:0*
T0*(
_output_shapes
:         а2
flatten/Reshapeб
dense/MatMul/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype02
dense/MatMul/ReadVariableOpШ
dense/MatMulMatMulflatten/Reshape:output:0#dense/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2
dense/MatMulЯ
dense/BiasAdd/ReadVariableOpReadVariableOp%dense_biasadd_readvariableop_resource*
_output_shapes	
:ю*
dtype02
dense/BiasAdd/ReadVariableOpЪ
dense/BiasAddBiasAdddense/MatMul:product:0$dense/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         ю2
dense/BiasAddk

dense/ReluReludense/BiasAdd:output:0*
T0*(
_output_shapes
:         ю2

dense/Reluз
dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype02
dense_1/MatMul/ReadVariableOpЮ
dense_1/MatMulMatMuldense/Relu:activations:0%dense_1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2
dense_1/MatMulе
dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes	
:м*
dtype02 
dense_1/BiasAdd/ReadVariableOpв
dense_1/BiasAddBiasAdddense_1/MatMul:product:0&dense_1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:         м2
dense_1/BiasAddq
dense_1/ReluReludense_1/BiasAdd:output:0*
T0*(
_output_shapes
:         м2
dense_1/Relus
dropout/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n█╢?2
dropout/dropout/Constа
dropout/dropout/MulMuldense_1/Relu:activations:0dropout/dropout/Const:output:0*
T0*(
_output_shapes
:         м2
dropout/dropout/Mulx
dropout/dropout/ShapeShapedense_1/Relu:activations:0*
T0*
_output_shapes
:2
dropout/dropout/Shape═
,dropout/dropout/random_uniform/RandomUniformRandomUniformdropout/dropout/Shape:output:0*
T0*(
_output_shapes
:         м*
dtype02.
,dropout/dropout/random_uniform/RandomUniformЕ
dropout/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЪЩЩ>2 
dropout/dropout/GreaterEqual/y▀
dropout/dropout/GreaterEqualGreaterEqual5dropout/dropout/random_uniform/RandomUniform:output:0'dropout/dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:         м2
dropout/dropout/GreaterEqualШ
dropout/dropout/CastCast dropout/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:         м2
dropout/dropout/CastЫ
dropout/dropout/Mul_1Muldropout/dropout/Mul:z:0dropout/dropout/Cast:y:0*
T0*(
_output_shapes
:         м2
dropout/dropout/Mul_1ж
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes
:	мd*
dtype02
dense_2/MatMul/ReadVariableOpЮ
dense_2/MatMulMatMuldropout/dropout/Mul_1:z:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2
dense_2/MatMulд
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02 
dense_2/BiasAdd/ReadVariableOpб
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2
dense_2/BiasAddp
dense_2/ReluReludense_2/BiasAdd:output:0*
T0*'
_output_shapes
:         d2
dense_2/Reluе
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

:d2*
dtype02
dense_3/MatMul/ReadVariableOpЯ
dense_3/MatMulMatMuldense_2/Relu:activations:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22
dense_3/MatMulд
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype02 
dense_3/BiasAdd/ReadVariableOpб
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         22
dense_3/BiasAddp
dense_3/ReluReludense_3/BiasAdd:output:0*
T0*'
_output_shapes
:         22
dense_3/Reluw
dropout_1/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n█╢?2
dropout_1/dropout/Constе
dropout_1/dropout/MulMuldense_3/Relu:activations:0 dropout_1/dropout/Const:output:0*
T0*'
_output_shapes
:         22
dropout_1/dropout/Mul|
dropout_1/dropout/ShapeShapedense_3/Relu:activations:0*
T0*
_output_shapes
:2
dropout_1/dropout/Shape╥
.dropout_1/dropout/random_uniform/RandomUniformRandomUniform dropout_1/dropout/Shape:output:0*
T0*'
_output_shapes
:         2*
dtype020
.dropout_1/dropout/random_uniform/RandomUniformЙ
 dropout_1/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *ЪЩЩ>2"
 dropout_1/dropout/GreaterEqual/yц
dropout_1/dropout/GreaterEqualGreaterEqual7dropout_1/dropout/random_uniform/RandomUniform:output:0)dropout_1/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:         22 
dropout_1/dropout/GreaterEqualЭ
dropout_1/dropout/CastCast"dropout_1/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:         22
dropout_1/dropout/Castв
dropout_1/dropout/Mul_1Muldropout_1/dropout/Mul:z:0dropout_1/dropout/Cast:y:0*
T0*'
_output_shapes
:         22
dropout_1/dropout/Mul_1в
output/MatMul/ReadVariableOpReadVariableOp%output_matmul_readvariableop_resource*
_output_shapes

:2*
dtype02
output/MatMul/ReadVariableOpЭ
output/MatMulMatMuldropout_1/dropout/Mul_1:z:0$output/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
output/MatMulб
output/BiasAdd/ReadVariableOpReadVariableOp&output_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
output/BiasAdd/ReadVariableOpЭ
output/BiasAddBiasAddoutput/MatMul:product:0%output/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
output/BiasAddv
output/SoftmaxSoftmaxoutput/BiasAdd:output:0*
T0*'
_output_shapes
:         2
output/Softmax╟
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource* 
_output_shapes
:
аю*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOpп
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
аю2!
dense/kernel/Regularizer/SquareС
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const▓
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/SumЕ
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2 
dense/kernel/Regularizer/mul/x┤
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul═
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource* 
_output_shapes
:
юм*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp╡
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0* 
_output_shapes
:
юм2#
!dense_1/kernel/Regularizer/SquareХ
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const║
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/SumЙ
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_1/kernel/Regularizer/mul/x╝
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul╠
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul╦
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

:d2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp│
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d22#
!dense_3/kernel/Regularizer/SquareХ
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const║
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/SumЙ
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_3/kernel/Regularizer/mul/x╝
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/muls
IdentityIdentityoutput/Softmax:softmax:0^NoOp*
T0*'
_output_shapes
:         2

Identity│

NoOpNoOp^conv2d/BiasAdd/ReadVariableOp^conv2d/Conv2D/ReadVariableOp ^conv2d_1/BiasAdd/ReadVariableOp^conv2d_1/Conv2D/ReadVariableOp ^conv2d_2/BiasAdd/ReadVariableOp^conv2d_2/Conv2D/ReadVariableOp ^conv2d_3/BiasAdd/ReadVariableOp^conv2d_3/Conv2D/ReadVariableOp ^conv2d_4/BiasAdd/ReadVariableOp^conv2d_4/Conv2D/ReadVariableOp^dense/BiasAdd/ReadVariableOp^dense/MatMul/ReadVariableOp/^dense/kernel/Regularizer/Square/ReadVariableOp^dense_1/BiasAdd/ReadVariableOp^dense_1/MatMul/ReadVariableOp1^dense_1/kernel/Regularizer/Square/ReadVariableOp^dense_2/BiasAdd/ReadVariableOp^dense_2/MatMul/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp^output/BiasAdd/ReadVariableOp^output/MatMul/ReadVariableOp@^sequential/random_flip/stateful_uniform_full_int/RngReadAndSkipg^sequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgn^sequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter;^sequential/random_rotation/stateful_uniform/RngReadAndSkip7^sequential/random_zoom/stateful_uniform/RngReadAndSkip*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*\
_input_shapesK
I:         &: : : : : : : : : : : : : : : : : : : : : : : 2>
conv2d/BiasAdd/ReadVariableOpconv2d/BiasAdd/ReadVariableOp2<
conv2d/Conv2D/ReadVariableOpconv2d/Conv2D/ReadVariableOp2B
conv2d_1/BiasAdd/ReadVariableOpconv2d_1/BiasAdd/ReadVariableOp2@
conv2d_1/Conv2D/ReadVariableOpconv2d_1/Conv2D/ReadVariableOp2B
conv2d_2/BiasAdd/ReadVariableOpconv2d_2/BiasAdd/ReadVariableOp2@
conv2d_2/Conv2D/ReadVariableOpconv2d_2/Conv2D/ReadVariableOp2B
conv2d_3/BiasAdd/ReadVariableOpconv2d_3/BiasAdd/ReadVariableOp2@
conv2d_3/Conv2D/ReadVariableOpconv2d_3/Conv2D/ReadVariableOp2B
conv2d_4/BiasAdd/ReadVariableOpconv2d_4/BiasAdd/ReadVariableOp2@
conv2d_4/Conv2D/ReadVariableOpconv2d_4/Conv2D/ReadVariableOp2<
dense/BiasAdd/ReadVariableOpdense/BiasAdd/ReadVariableOp2:
dense/MatMul/ReadVariableOpdense/MatMul/ReadVariableOp2`
.dense/kernel/Regularizer/Square/ReadVariableOp.dense/kernel/Regularizer/Square/ReadVariableOp2@
dense_1/BiasAdd/ReadVariableOpdense_1/BiasAdd/ReadVariableOp2>
dense_1/MatMul/ReadVariableOpdense_1/MatMul/ReadVariableOp2d
0dense_1/kernel/Regularizer/Square/ReadVariableOp0dense_1/kernel/Regularizer/Square/ReadVariableOp2@
dense_2/BiasAdd/ReadVariableOpdense_2/BiasAdd/ReadVariableOp2>
dense_2/MatMul/ReadVariableOpdense_2/MatMul/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2>
output/BiasAdd/ReadVariableOpoutput/BiasAdd/ReadVariableOp2<
output/MatMul/ReadVariableOpoutput/MatMul/ReadVariableOp2В
?sequential/random_flip/stateful_uniform_full_int/RngReadAndSkip?sequential/random_flip/stateful_uniform_full_int/RngReadAndSkip2╨
fsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlgfsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetAlg2▐
msequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCountermsequential/random_flip/stateless_random_flip_left_right/stateless_random_uniform/StatelessRandomGetKeyCounter2x
:sequential/random_rotation/stateful_uniform/RngReadAndSkip:sequential/random_rotation/stateful_uniform/RngReadAndSkip2p
6sequential/random_zoom/stateful_uniform/RngReadAndSkip6sequential/random_zoom/stateful_uniform/RngReadAndSkip:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
╜
f
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_78324

inputs
identityТ
MaxPoolMaxPoolinputs*/
_output_shapes
:         	*
ksize
*
paddingVALID*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:         	2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         :W S
/
_output_shapes
:         
 
_user_specified_nameinputs
щ
№
C__inference_conv2d_2_layer_call_and_return_conditional_losses_78344

inputs8
conv2d_readvariableop_resource: -
biasadd_readvariableop_resource: 
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: *
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 *
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         	 2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         	 2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         	 2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         	: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:         	
 
_user_specified_nameinputs
к
f
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_78419

inputs
identityн
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4                                    *
ksize
*
paddingVALID*
strides
2	
MaxPoolЗ
IdentityIdentityMaxPool:output:0*
T0*J
_output_shapes8
6:4                                    2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4                                    :r n
J
_output_shapes8
6:4                                    
 
_user_specified_nameinputs
┌
┴
,__inference_sequential_1_layer_call_fn_76638
sequential_input!
unknown:
	unknown_0:#
	unknown_1:
	unknown_2:#
	unknown_3: 
	unknown_4: #
	unknown_5: @
	unknown_6:@#
	unknown_7:@P
	unknown_8:P
	unknown_9:
аю

unknown_10:	ю

unknown_11:
юм

unknown_12:	м

unknown_13:	мd

unknown_14:d

unknown_15:d2

unknown_16:2

unknown_17:2

unknown_18:
identityИвStatefulPartitionedCallЇ
StatefulPartitionedCallStatefulPartitionedCallsequential_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *6
_read_only_resource_inputs
	
*-
config_proto

CPU

GPU 2J 8В *P
fKRI
G__inference_sequential_1_layer_call_and_return_conditional_losses_765952
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*V
_input_shapesE
C:         &: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:a ]
/
_output_shapes
:         &
*
_user_specified_namesequential_input
┤С
┐
F__inference_random_zoom_layer_call_and_return_conditional_losses_79023

inputs6
(stateful_uniform_rngreadandskip_resource:	
identityИвstateful_uniform/RngReadAndSkipD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliceБ
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:
¤        2
strided_slice_1/stackЕ
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ь
strided_slice_1StridedSliceShape:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1^
CastCaststrided_slice_1:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
CastБ
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice_2/stackЕ
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ь
strided_slice_2StridedSliceShape:output:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_2b
Cast_1Caststrided_slice_2:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
Cast_1v
stateful_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2
stateful_uniform/shape/1б
stateful_uniform/shapePackstrided_slice:output:0!stateful_uniform/shape/1:output:0*
N*
T0*
_output_shapes
:2
stateful_uniform/shapeq
stateful_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *fff?2
stateful_uniform/minq
stateful_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *═╠М?2
stateful_uniform/maxz
stateful_uniform/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
stateful_uniform/ConstЩ
stateful_uniform/ProdProdstateful_uniform/shape:output:0stateful_uniform/Const:output:0*
T0*
_output_shapes
: 2
stateful_uniform/Prodt
stateful_uniform/Cast/xConst*
_output_shapes
: *
dtype0*
value	B :2
stateful_uniform/Cast/xК
stateful_uniform/Cast_1Caststateful_uniform/Prod:output:0*

DstT0*

SrcT0*
_output_shapes
: 2
stateful_uniform/Cast_1┘
stateful_uniform/RngReadAndSkipRngReadAndSkip(stateful_uniform_rngreadandskip_resource stateful_uniform/Cast/x:output:0stateful_uniform/Cast_1:y:0*
_output_shapes
:2!
stateful_uniform/RngReadAndSkipЦ
$stateful_uniform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2&
$stateful_uniform/strided_slice/stackЪ
&stateful_uniform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice/stack_1Ъ
&stateful_uniform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice/stack_2╬
stateful_uniform/strided_sliceStridedSlice'stateful_uniform/RngReadAndSkip:value:0-stateful_uniform/strided_slice/stack:output:0/stateful_uniform/strided_slice/stack_1:output:0/stateful_uniform/strided_slice/stack_2:output:0*
Index0*
T0	*
_output_shapes
:*

begin_mask2 
stateful_uniform/strided_sliceЩ
stateful_uniform/BitcastBitcast'stateful_uniform/strided_slice:output:0*
T0	*
_output_shapes
:*

type02
stateful_uniform/BitcastЪ
&stateful_uniform/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB:2(
&stateful_uniform/strided_slice_1/stackЮ
(stateful_uniform/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(stateful_uniform/strided_slice_1/stack_1Ю
(stateful_uniform/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(stateful_uniform/strided_slice_1/stack_2╞
 stateful_uniform/strided_slice_1StridedSlice'stateful_uniform/RngReadAndSkip:value:0/stateful_uniform/strided_slice_1/stack:output:01stateful_uniform/strided_slice_1/stack_1:output:01stateful_uniform/strided_slice_1/stack_2:output:0*
Index0*
T0	*
_output_shapes
:2"
 stateful_uniform/strided_slice_1Я
stateful_uniform/Bitcast_1Bitcast)stateful_uniform/strided_slice_1:output:0*
T0	*
_output_shapes
:*

type02
stateful_uniform/Bitcast_1а
-stateful_uniform/StatelessRandomUniformV2/algConst*
_output_shapes
: *
dtype0*
value	B :2/
-stateful_uniform/StatelessRandomUniformV2/alg╝
)stateful_uniform/StatelessRandomUniformV2StatelessRandomUniformV2stateful_uniform/shape:output:0#stateful_uniform/Bitcast_1:output:0!stateful_uniform/Bitcast:output:06stateful_uniform/StatelessRandomUniformV2/alg:output:0*'
_output_shapes
:         2+
)stateful_uniform/StatelessRandomUniformV2Т
stateful_uniform/subSubstateful_uniform/max:output:0stateful_uniform/min:output:0*
T0*
_output_shapes
: 2
stateful_uniform/sub│
stateful_uniform/mulMul2stateful_uniform/StatelessRandomUniformV2:output:0stateful_uniform/sub:z:0*
T0*'
_output_shapes
:         2
stateful_uniform/mulШ
stateful_uniformAddV2stateful_uniform/mul:z:0stateful_uniform/min:output:0*
T0*'
_output_shapes
:         2
stateful_uniform\
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concat/axisЩ
concatConcatV2stateful_uniform:z:0stateful_uniform:z:0concat/axis:output:0*
N*
T0*'
_output_shapes
:         2
concate
zoom_matrix/ShapeShapeconcat:output:0*
T0*
_output_shapes
:2
zoom_matrix/ShapeМ
zoom_matrix/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2!
zoom_matrix/strided_slice/stackР
!zoom_matrix/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2#
!zoom_matrix/strided_slice/stack_1Р
!zoom_matrix/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2#
!zoom_matrix/strided_slice/stack_2к
zoom_matrix/strided_sliceStridedSlicezoom_matrix/Shape:output:0(zoom_matrix/strided_slice/stack:output:0*zoom_matrix/strided_slice/stack_1:output:0*zoom_matrix/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
zoom_matrix/strided_slicek
zoom_matrix/sub/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
zoom_matrix/sub/yr
zoom_matrix/subSub
Cast_1:y:0zoom_matrix/sub/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/subs
zoom_matrix/truediv/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
zoom_matrix/truediv/yЛ
zoom_matrix/truedivRealDivzoom_matrix/sub:z:0zoom_matrix/truediv/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/truedivЫ
!zoom_matrix/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2#
!zoom_matrix/strided_slice_1/stackЯ
#zoom_matrix/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2%
#zoom_matrix/strided_slice_1/stack_1Я
#zoom_matrix/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2%
#zoom_matrix/strided_slice_1/stack_2ё
zoom_matrix/strided_slice_1StridedSliceconcat:output:0*zoom_matrix/strided_slice_1/stack:output:0,zoom_matrix/strided_slice_1/stack_1:output:0,zoom_matrix/strided_slice_1/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2
zoom_matrix/strided_slice_1o
zoom_matrix/sub_1/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
zoom_matrix/sub_1/xг
zoom_matrix/sub_1Subzoom_matrix/sub_1/x:output:0$zoom_matrix/strided_slice_1:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/sub_1Л
zoom_matrix/mulMulzoom_matrix/truediv:z:0zoom_matrix/sub_1:z:0*
T0*'
_output_shapes
:         2
zoom_matrix/mulo
zoom_matrix/sub_2/yConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
zoom_matrix/sub_2/yv
zoom_matrix/sub_2SubCast:y:0zoom_matrix/sub_2/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/sub_2w
zoom_matrix/truediv_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
zoom_matrix/truediv_1/yУ
zoom_matrix/truediv_1RealDivzoom_matrix/sub_2:z:0 zoom_matrix/truediv_1/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/truediv_1Ы
!zoom_matrix/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2#
!zoom_matrix/strided_slice_2/stackЯ
#zoom_matrix/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2%
#zoom_matrix/strided_slice_2/stack_1Я
#zoom_matrix/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2%
#zoom_matrix/strided_slice_2/stack_2ё
zoom_matrix/strided_slice_2StridedSliceconcat:output:0*zoom_matrix/strided_slice_2/stack:output:0,zoom_matrix/strided_slice_2/stack_1:output:0,zoom_matrix/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2
zoom_matrix/strided_slice_2o
zoom_matrix/sub_3/xConst*
_output_shapes
: *
dtype0*
valueB
 *  А?2
zoom_matrix/sub_3/xг
zoom_matrix/sub_3Subzoom_matrix/sub_3/x:output:0$zoom_matrix/strided_slice_2:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/sub_3С
zoom_matrix/mul_1Mulzoom_matrix/truediv_1:z:0zoom_matrix/sub_3:z:0*
T0*'
_output_shapes
:         2
zoom_matrix/mul_1Ы
!zoom_matrix/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2#
!zoom_matrix/strided_slice_3/stackЯ
#zoom_matrix/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2%
#zoom_matrix/strided_slice_3/stack_1Я
#zoom_matrix/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2%
#zoom_matrix/strided_slice_3/stack_2ё
zoom_matrix/strided_slice_3StridedSliceconcat:output:0*zoom_matrix/strided_slice_3/stack:output:0,zoom_matrix/strided_slice_3/stack_1:output:0,zoom_matrix/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2
zoom_matrix/strided_slice_3t
zoom_matrix/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros/mul/yЬ
zoom_matrix/zeros/mulMul"zoom_matrix/strided_slice:output:0 zoom_matrix/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros/mulw
zoom_matrix/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zoom_matrix/zeros/Less/yЧ
zoom_matrix/zeros/LessLesszoom_matrix/zeros/mul:z:0!zoom_matrix/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros/Lessz
zoom_matrix/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros/packed/1│
zoom_matrix/zeros/packedPack"zoom_matrix/strided_slice:output:0#zoom_matrix/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zoom_matrix/zeros/packedw
zoom_matrix/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zoom_matrix/zeros/Constе
zoom_matrix/zerosFill!zoom_matrix/zeros/packed:output:0 zoom_matrix/zeros/Const:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/zerosx
zoom_matrix/zeros_1/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros_1/mul/yв
zoom_matrix/zeros_1/mulMul"zoom_matrix/strided_slice:output:0"zoom_matrix/zeros_1/mul/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros_1/mul{
zoom_matrix/zeros_1/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zoom_matrix/zeros_1/Less/yЯ
zoom_matrix/zeros_1/LessLesszoom_matrix/zeros_1/mul:z:0#zoom_matrix/zeros_1/Less/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros_1/Less~
zoom_matrix/zeros_1/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros_1/packed/1╣
zoom_matrix/zeros_1/packedPack"zoom_matrix/strided_slice:output:0%zoom_matrix/zeros_1/packed/1:output:0*
N*
T0*
_output_shapes
:2
zoom_matrix/zeros_1/packed{
zoom_matrix/zeros_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zoom_matrix/zeros_1/Constн
zoom_matrix/zeros_1Fill#zoom_matrix/zeros_1/packed:output:0"zoom_matrix/zeros_1/Const:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/zeros_1Ы
!zoom_matrix/strided_slice_4/stackConst*
_output_shapes
:*
dtype0*!
valueB"           2#
!zoom_matrix/strided_slice_4/stackЯ
#zoom_matrix/strided_slice_4/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"           2%
#zoom_matrix/strided_slice_4/stack_1Я
#zoom_matrix/strided_slice_4/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2%
#zoom_matrix/strided_slice_4/stack_2ё
zoom_matrix/strided_slice_4StridedSliceconcat:output:0*zoom_matrix/strided_slice_4/stack:output:0,zoom_matrix/strided_slice_4/stack_1:output:0,zoom_matrix/strided_slice_4/stack_2:output:0*
Index0*
T0*'
_output_shapes
:         *

begin_mask*
end_mask*
new_axis_mask*
shrink_axis_mask2
zoom_matrix/strided_slice_4x
zoom_matrix/zeros_2/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros_2/mul/yв
zoom_matrix/zeros_2/mulMul"zoom_matrix/strided_slice:output:0"zoom_matrix/zeros_2/mul/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros_2/mul{
zoom_matrix/zeros_2/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zoom_matrix/zeros_2/Less/yЯ
zoom_matrix/zeros_2/LessLesszoom_matrix/zeros_2/mul:z:0#zoom_matrix/zeros_2/Less/y:output:0*
T0*
_output_shapes
: 2
zoom_matrix/zeros_2/Less~
zoom_matrix/zeros_2/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/zeros_2/packed/1╣
zoom_matrix/zeros_2/packedPack"zoom_matrix/strided_slice:output:0%zoom_matrix/zeros_2/packed/1:output:0*
N*
T0*
_output_shapes
:2
zoom_matrix/zeros_2/packed{
zoom_matrix/zeros_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zoom_matrix/zeros_2/Constн
zoom_matrix/zeros_2Fill#zoom_matrix/zeros_2/packed:output:0"zoom_matrix/zeros_2/Const:output:0*
T0*'
_output_shapes
:         2
zoom_matrix/zeros_2t
zoom_matrix/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
zoom_matrix/concat/axisс
zoom_matrix/concatConcatV2$zoom_matrix/strided_slice_3:output:0zoom_matrix/zeros:output:0zoom_matrix/mul:z:0zoom_matrix/zeros_1:output:0$zoom_matrix/strided_slice_4:output:0zoom_matrix/mul_1:z:0zoom_matrix/zeros_2:output:0 zoom_matrix/concat/axis:output:0*
N*
T0*'
_output_shapes
:         2
zoom_matrix/concatX
transform/ShapeShapeinputs*
T0*
_output_shapes
:2
transform/ShapeИ
transform/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2
transform/strided_slice/stackМ
transform/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2!
transform/strided_slice/stack_1М
transform/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2!
transform/strided_slice/stack_2К
transform/strided_sliceStridedSlicetransform/Shape:output:0&transform/strided_slice/stack:output:0(transform/strided_slice/stack_1:output:0(transform/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2
transform/strided_sliceq
transform/fill_valueConst*
_output_shapes
: *
dtype0*
valueB
 *    2
transform/fill_value├
$transform/ImageProjectiveTransformV3ImageProjectiveTransformV3inputszoom_matrix/concat:output:0 transform/strided_slice:output:0transform/fill_value:output:0*/
_output_shapes
:         &*
dtype0*
	fill_mode	REFLECT*
interpolation
BILINEAR2&
$transform/ImageProjectiveTransformV3Ь
IdentityIdentity9transform/ImageProjectiveTransformV3:transformed_images:0^NoOp*
T0*/
_output_shapes
:         &2

Identityp
NoOpNoOp ^stateful_uniform/RngReadAndSkip*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 2B
stateful_uniform/RngReadAndSkipstateful_uniform/RngReadAndSkip:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
ш
K
/__inference_random_rotation_layer_call_fn_78760

inputs
identity╨
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *S
fNRL
J__inference_random_rotation_layer_call_and_return_conditional_losses_758192
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
ч
·
A__inference_conv2d_layer_call_and_return_conditional_losses_78264

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identityИвBiasAdd/ReadVariableOpвConv2D/ReadVariableOpХ
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOpг
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &*
paddingSAME*
strides
2
Conv2DМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpИ
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         &2	
BiasAdd`
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:         &2
Reluu
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:         &2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         &: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
ё
b
D__inference_dropout_1_layer_call_and_return_conditional_losses_78605

inputs

identity_1Z
IdentityIdentityinputs*
T0*'
_output_shapes
:         22

Identityi

Identity_1IdentityIdentity:output:0*
T0*'
_output_shapes
:         22

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*&
_input_shapes
:         2:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
р
G
+__inference_rescaling_1_layer_call_fn_78236

inputs
identity╠
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_rescaling_1_layer_call_and_return_conditional_losses_763282
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:         &2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:         &:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
я
з
B__inference_dense_2_layer_call_and_return_conditional_losses_76517

inputs1
matmul_readvariableop_resource:	мd-
biasadd_readvariableop_resource:d
identityИвBiasAdd/ReadVariableOpвMatMul/ReadVariableOpв0dense_2/kernel/Regularizer/Square/ReadVariableOpО
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	мd*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2
MatMulМ
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
BiasAdd/ReadVariableOpБ
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:         d2	
BiasAddX
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:         d2
Relu─
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	мd*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp┤
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	мd2#
!dense_2/kernel/Regularizer/SquareХ
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const║
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/SumЙ
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *oГ:2"
 dense_2/kernel/Regularizer/mul/x╝
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mulm
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:         d2

Identity▓
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:         м: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:         м
 
_user_specified_nameinputs
ц
{
+__inference_random_zoom_layer_call_fn_78905

inputs
unknown:	
identityИвStatefulPartitionedCallю
StatefulPartitionedCallStatefulPartitionedCallinputsunknown*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8В *O
fJRH
F__inference_random_zoom_layer_call_and_return_conditional_losses_759532
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         &2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*0
_input_shapes
:         &: 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs
ы
У
&__inference_output_layer_call_fn_78626

inputs
unknown:2
	unknown_0:
identityИвStatefulPartitionedCallё
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:         *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_output_layer_call_and_return_conditional_losses_765642
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:         2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:         2: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:         2
 
_user_specified_nameinputs
Ф
Ы
&__inference_conv2d_layer_call_fn_78253

inputs!
unknown:
	unknown_0:
identityИвStatefulPartitionedCall∙
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:         &*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8В *J
fERC
A__inference_conv2d_layer_call_and_return_conditional_losses_763412
StatefulPartitionedCallГ
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:         &2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:         &: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:         &
 
_user_specified_nameinputs"иL
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*├
serving_defaultп
U
sequential_inputA
"serving_default_sequential_input:0         &:
output0
StatefulPartitionedCall:0         tensorflow/serving/predict:ШЖ
э
layer-0
layer-1
layer_with_weights-0
layer-2
layer-3
layer_with_weights-1
layer-4
layer-5
layer_with_weights-2
layer-6
layer-7
	layer_with_weights-3
	layer-8

layer_with_weights-4

layer-9
layer-10
layer-11
layer_with_weights-5
layer-12
layer_with_weights-6
layer-13
layer-14
layer_with_weights-7
layer-15
layer_with_weights-8
layer-16
layer-17
layer_with_weights-9
layer-18
	optimizer
	variables
regularization_losses
trainable_variables
	keras_api

signatures
║__call__
╗_default_save_signature
+╝&call_and_return_all_conditional_losses"
_tf_keras_sequential
╙
layer-0
layer-1
layer-2
	variables
regularization_losses
trainable_variables
 	keras_api
╜__call__
+╛&call_and_return_all_conditional_losses"
_tf_keras_sequential
з
!	variables
"regularization_losses
#trainable_variables
$	keras_api
┐__call__
+└&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

%kernel
&bias
'	variables
(regularization_losses
)trainable_variables
*	keras_api
┴__call__
+┬&call_and_return_all_conditional_losses"
_tf_keras_layer
з
+	variables
,regularization_losses
-trainable_variables
.	keras_api
├__call__
+─&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

/kernel
0bias
1	variables
2regularization_losses
3trainable_variables
4	keras_api
┼__call__
+╞&call_and_return_all_conditional_losses"
_tf_keras_layer
з
5	variables
6regularization_losses
7trainable_variables
8	keras_api
╟__call__
+╚&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

9kernel
:bias
;	variables
<regularization_losses
=trainable_variables
>	keras_api
╔__call__
+╩&call_and_return_all_conditional_losses"
_tf_keras_layer
з
?	variables
@regularization_losses
Atrainable_variables
B	keras_api
╦__call__
+╠&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

Ckernel
Dbias
E	variables
Fregularization_losses
Gtrainable_variables
H	keras_api
═__call__
+╬&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

Ikernel
Jbias
K	variables
Lregularization_losses
Mtrainable_variables
N	keras_api
╧__call__
+╨&call_and_return_all_conditional_losses"
_tf_keras_layer
з
O	variables
Pregularization_losses
Qtrainable_variables
R	keras_api
╤__call__
+╥&call_and_return_all_conditional_losses"
_tf_keras_layer
з
S	variables
Tregularization_losses
Utrainable_variables
V	keras_api
╙__call__
+╘&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

Wkernel
Xbias
Y	variables
Zregularization_losses
[trainable_variables
\	keras_api
╒__call__
+╓&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

]kernel
^bias
_	variables
`regularization_losses
atrainable_variables
b	keras_api
╫__call__
+╪&call_and_return_all_conditional_losses"
_tf_keras_layer
з
c	variables
dregularization_losses
etrainable_variables
f	keras_api
┘__call__
+┌&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

gkernel
hbias
i	variables
jregularization_losses
ktrainable_variables
l	keras_api
█__call__
+▄&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

mkernel
nbias
o	variables
pregularization_losses
qtrainable_variables
r	keras_api
▌__call__
+▐&call_and_return_all_conditional_losses"
_tf_keras_layer
з
s	variables
tregularization_losses
utrainable_variables
v	keras_api
▀__call__
+р&call_and_return_all_conditional_losses"
_tf_keras_layer
╜

wkernel
xbias
y	variables
zregularization_losses
{trainable_variables
|	keras_api
с__call__
+т&call_and_return_all_conditional_losses"
_tf_keras_layer
х
}iter

~beta_1

beta_2

Аdecay
Бlearning_rate%mТ&mУ/mФ0mХ9mЦ:mЧCmШDmЩImЪJmЫWmЬXmЭ]mЮ^mЯgmаhmбmmвnmгwmдxmе%vж&vз/vи0vй9vк:vлCvмDvнIvоJvпWv░Xv▒]v▓^v│gv┤hv╡mv╢nv╖wv╕xv╣"
	optimizer
╢
%0
&1
/2
03
94
:5
C6
D7
I8
J9
W10
X11
]12
^13
g14
h15
m16
n17
w18
x19"
trackable_list_wrapper
@
у0
ф1
х2
ц3"
trackable_list_wrapper
╢
%0
&1
/2
03
94
:5
C6
D7
I8
J9
W10
X11
]12
^13
g14
h15
m16
n17
w18
x19"
trackable_list_wrapper
╙
Вmetrics
Гlayers
	variables
Дlayer_metrics
Еnon_trainable_variables
regularization_losses
 Жlayer_regularization_losses
trainable_variables
║__call__
╗_default_save_signature
+╝&call_and_return_all_conditional_losses
'╝"call_and_return_conditional_losses"
_generic_user_object
-
чserving_default"
signature_map
╢
	З_rng
И	variables
Йregularization_losses
Кtrainable_variables
Л	keras_api
ш__call__
+щ&call_and_return_all_conditional_losses"
_tf_keras_layer
╢
	М_rng
Н	variables
Оregularization_losses
Пtrainable_variables
Р	keras_api
ъ__call__
+ы&call_and_return_all_conditional_losses"
_tf_keras_layer
╢
	С_rng
Т	variables
Уregularization_losses
Фtrainable_variables
Х	keras_api
ь__call__
+э&call_and_return_all_conditional_losses"
_tf_keras_layer
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
Цmetrics
Чlayers
	variables
Шlayer_metrics
Щnon_trainable_variables
regularization_losses
 Ъlayer_regularization_losses
trainable_variables
╜__call__
+╛&call_and_return_all_conditional_losses
'╛"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
Ыmetrics
Ьlayers
Эlayer_metrics
!	variables
Юnon_trainable_variables
"regularization_losses
 Яlayer_regularization_losses
#trainable_variables
┐__call__
+└&call_and_return_all_conditional_losses
'└"call_and_return_conditional_losses"
_generic_user_object
':%2conv2d/kernel
:2conv2d/bias
.
%0
&1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
╡
аmetrics
бlayers
вlayer_metrics
'	variables
гnon_trainable_variables
(regularization_losses
 дlayer_regularization_losses
)trainable_variables
┴__call__
+┬&call_and_return_all_conditional_losses
'┬"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
еmetrics
жlayers
зlayer_metrics
+	variables
иnon_trainable_variables
,regularization_losses
 йlayer_regularization_losses
-trainable_variables
├__call__
+─&call_and_return_all_conditional_losses
'─"call_and_return_conditional_losses"
_generic_user_object
):'2conv2d_1/kernel
:2conv2d_1/bias
.
/0
01"
trackable_list_wrapper
 "
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
╡
кmetrics
лlayers
мlayer_metrics
1	variables
нnon_trainable_variables
2regularization_losses
 оlayer_regularization_losses
3trainable_variables
┼__call__
+╞&call_and_return_all_conditional_losses
'╞"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
пmetrics
░layers
▒layer_metrics
5	variables
▓non_trainable_variables
6regularization_losses
 │layer_regularization_losses
7trainable_variables
╟__call__
+╚&call_and_return_all_conditional_losses
'╚"call_and_return_conditional_losses"
_generic_user_object
):' 2conv2d_2/kernel
: 2conv2d_2/bias
.
90
:1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
90
:1"
trackable_list_wrapper
╡
┤metrics
╡layers
╢layer_metrics
;	variables
╖non_trainable_variables
<regularization_losses
 ╕layer_regularization_losses
=trainable_variables
╔__call__
+╩&call_and_return_all_conditional_losses
'╩"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
╣metrics
║layers
╗layer_metrics
?	variables
╝non_trainable_variables
@regularization_losses
 ╜layer_regularization_losses
Atrainable_variables
╦__call__
+╠&call_and_return_all_conditional_losses
'╠"call_and_return_conditional_losses"
_generic_user_object
):' @2conv2d_3/kernel
:@2conv2d_3/bias
.
C0
D1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
C0
D1"
trackable_list_wrapper
╡
╛metrics
┐layers
└layer_metrics
E	variables
┴non_trainable_variables
Fregularization_losses
 ┬layer_regularization_losses
Gtrainable_variables
═__call__
+╬&call_and_return_all_conditional_losses
'╬"call_and_return_conditional_losses"
_generic_user_object
):'@P2conv2d_4/kernel
:P2conv2d_4/bias
.
I0
J1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
I0
J1"
trackable_list_wrapper
╡
├metrics
─layers
┼layer_metrics
K	variables
╞non_trainable_variables
Lregularization_losses
 ╟layer_regularization_losses
Mtrainable_variables
╧__call__
+╨&call_and_return_all_conditional_losses
'╨"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
╚metrics
╔layers
╩layer_metrics
O	variables
╦non_trainable_variables
Pregularization_losses
 ╠layer_regularization_losses
Qtrainable_variables
╤__call__
+╥&call_and_return_all_conditional_losses
'╥"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
═metrics
╬layers
╧layer_metrics
S	variables
╨non_trainable_variables
Tregularization_losses
 ╤layer_regularization_losses
Utrainable_variables
╙__call__
+╘&call_and_return_all_conditional_losses
'╘"call_and_return_conditional_losses"
_generic_user_object
 :
аю2dense/kernel
:ю2
dense/bias
.
W0
X1"
trackable_list_wrapper
(
у0"
trackable_list_wrapper
.
W0
X1"
trackable_list_wrapper
╡
╥metrics
╙layers
╘layer_metrics
Y	variables
╒non_trainable_variables
Zregularization_losses
 ╓layer_regularization_losses
[trainable_variables
╒__call__
+╓&call_and_return_all_conditional_losses
'╓"call_and_return_conditional_losses"
_generic_user_object
": 
юм2dense_1/kernel
:м2dense_1/bias
.
]0
^1"
trackable_list_wrapper
(
ф0"
trackable_list_wrapper
.
]0
^1"
trackable_list_wrapper
╡
╫metrics
╪layers
┘layer_metrics
_	variables
┌non_trainable_variables
`regularization_losses
 █layer_regularization_losses
atrainable_variables
╫__call__
+╪&call_and_return_all_conditional_losses
'╪"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
▄metrics
▌layers
▐layer_metrics
c	variables
▀non_trainable_variables
dregularization_losses
 рlayer_regularization_losses
etrainable_variables
┘__call__
+┌&call_and_return_all_conditional_losses
'┌"call_and_return_conditional_losses"
_generic_user_object
!:	мd2dense_2/kernel
:d2dense_2/bias
.
g0
h1"
trackable_list_wrapper
(
х0"
trackable_list_wrapper
.
g0
h1"
trackable_list_wrapper
╡
сmetrics
тlayers
уlayer_metrics
i	variables
фnon_trainable_variables
jregularization_losses
 хlayer_regularization_losses
ktrainable_variables
█__call__
+▄&call_and_return_all_conditional_losses
'▄"call_and_return_conditional_losses"
_generic_user_object
 :d22dense_3/kernel
:22dense_3/bias
.
m0
n1"
trackable_list_wrapper
(
ц0"
trackable_list_wrapper
.
m0
n1"
trackable_list_wrapper
╡
цmetrics
чlayers
шlayer_metrics
o	variables
щnon_trainable_variables
pregularization_losses
 ъlayer_regularization_losses
qtrainable_variables
▌__call__
+▐&call_and_return_all_conditional_losses
'▐"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╡
ыmetrics
ьlayers
эlayer_metrics
s	variables
юnon_trainable_variables
tregularization_losses
 яlayer_regularization_losses
utrainable_variables
▀__call__
+р&call_and_return_all_conditional_losses
'р"call_and_return_conditional_losses"
_generic_user_object
:22output/kernel
:2output/bias
.
w0
x1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
w0
x1"
trackable_list_wrapper
╡
Ёmetrics
ёlayers
Єlayer_metrics
y	variables
єnon_trainable_variables
zregularization_losses
 Їlayer_regularization_losses
{trainable_variables
с__call__
+т&call_and_return_all_conditional_losses
'т"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
0
ї0
Ў1"
trackable_list_wrapper
о
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
14
15
16
17
18"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
/
ў
_state_var"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╕
°metrics
∙layers
·layer_metrics
И	variables
√non_trainable_variables
Йregularization_losses
 №layer_regularization_losses
Кtrainable_variables
ш__call__
+щ&call_and_return_all_conditional_losses
'щ"call_and_return_conditional_losses"
_generic_user_object
/
¤
_state_var"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╕
■metrics
 layers
Аlayer_metrics
Н	variables
Бnon_trainable_variables
Оregularization_losses
 Вlayer_regularization_losses
Пtrainable_variables
ъ__call__
+ы&call_and_return_all_conditional_losses
'ы"call_and_return_conditional_losses"
_generic_user_object
/
Г
_state_var"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╕
Дmetrics
Еlayers
Жlayer_metrics
Т	variables
Зnon_trainable_variables
Уregularization_losses
 Иlayer_regularization_losses
Фtrainable_variables
ь__call__
+э&call_and_return_all_conditional_losses
'э"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
5
0
1
2"
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
(
у0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
(
ф0"
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
(
х0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
(
ц0"
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
R

Йtotal

Кcount
Л	variables
М	keras_api"
_tf_keras_metric
c

Нtotal

Оcount
П
_fn_kwargs
Р	variables
С	keras_api"
_tf_keras_metric
:	2Variable
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
:	2Variable
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
:	2Variable
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
:  (2total
:  (2count
0
Й0
К1"
trackable_list_wrapper
.
Л	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
0
Н0
О1"
trackable_list_wrapper
.
Р	variables"
_generic_user_object
,:*2Adam/conv2d/kernel/m
:2Adam/conv2d/bias/m
.:,2Adam/conv2d_1/kernel/m
 :2Adam/conv2d_1/bias/m
.:, 2Adam/conv2d_2/kernel/m
 : 2Adam/conv2d_2/bias/m
.:, @2Adam/conv2d_3/kernel/m
 :@2Adam/conv2d_3/bias/m
.:,@P2Adam/conv2d_4/kernel/m
 :P2Adam/conv2d_4/bias/m
%:#
аю2Adam/dense/kernel/m
:ю2Adam/dense/bias/m
':%
юм2Adam/dense_1/kernel/m
 :м2Adam/dense_1/bias/m
&:$	мd2Adam/dense_2/kernel/m
:d2Adam/dense_2/bias/m
%:#d22Adam/dense_3/kernel/m
:22Adam/dense_3/bias/m
$:"22Adam/output/kernel/m
:2Adam/output/bias/m
,:*2Adam/conv2d/kernel/v
:2Adam/conv2d/bias/v
.:,2Adam/conv2d_1/kernel/v
 :2Adam/conv2d_1/bias/v
.:, 2Adam/conv2d_2/kernel/v
 : 2Adam/conv2d_2/bias/v
.:, @2Adam/conv2d_3/kernel/v
 :@2Adam/conv2d_3/bias/v
.:,@P2Adam/conv2d_4/kernel/v
 :P2Adam/conv2d_4/bias/v
%:#
аю2Adam/dense/kernel/v
:ю2Adam/dense/bias/v
':%
юм2Adam/dense_1/kernel/v
 :м2Adam/dense_1/bias/v
&:$	мd2Adam/dense_2/kernel/v
:d2Adam/dense_2/bias/v
%:#d22Adam/dense_3/kernel/v
:22Adam/dense_3/bias/v
$:"22Adam/output/kernel/v
:2Adam/output/bias/v
■2√
,__inference_sequential_1_layer_call_fn_76638
,__inference_sequential_1_layer_call_fn_77358
,__inference_sequential_1_layer_call_fn_77409
,__inference_sequential_1_layer_call_fn_77056└
╖▓│
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╘B╤
 __inference__wrapped_model_75802sequential_input"Ш
С▓Н
FullArgSpec
argsЪ 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ъ2ч
G__inference_sequential_1_layer_call_and_return_conditional_losses_77519
G__inference_sequential_1_layer_call_and_return_conditional_losses_77925
G__inference_sequential_1_layer_call_and_return_conditional_losses_77143
G__inference_sequential_1_layer_call_and_return_conditional_losses_77236└
╖▓│
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
Ў2є
*__inference_sequential_layer_call_fn_75831
*__inference_sequential_layer_call_fn_77930
*__inference_sequential_layer_call_fn_77941
*__inference_sequential_layer_call_fn_76204└
╖▓│
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
т2▀
E__inference_sequential_layer_call_and_return_conditional_losses_77945
E__inference_sequential_layer_call_and_return_conditional_losses_78231
E__inference_sequential_layer_call_and_return_conditional_losses_76211
E__inference_sequential_layer_call_and_return_conditional_losses_76224└
╖▓│
FullArgSpec1
args)Ъ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsЪ
p 

 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╒2╥
+__inference_rescaling_1_layer_call_fn_78236в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
Ё2э
F__inference_rescaling_1_layer_call_and_return_conditional_losses_78244в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╨2═
&__inference_conv2d_layer_call_fn_78253в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_conv2d_layer_call_and_return_conditional_losses_78264в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
Ж2Г
-__inference_max_pooling2d_layer_call_fn_78269
-__inference_max_pooling2d_layer_call_fn_78274в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╝2╣
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_78279
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_78284в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╥2╧
(__inference_conv2d_1_layer_call_fn_78293в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
э2ъ
C__inference_conv2d_1_layer_call_and_return_conditional_losses_78304в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
К2З
/__inference_max_pooling2d_1_layer_call_fn_78309
/__inference_max_pooling2d_1_layer_call_fn_78314в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
└2╜
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_78319
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_78324в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╥2╧
(__inference_conv2d_2_layer_call_fn_78333в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
э2ъ
C__inference_conv2d_2_layer_call_and_return_conditional_losses_78344в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
К2З
/__inference_max_pooling2d_2_layer_call_fn_78349
/__inference_max_pooling2d_2_layer_call_fn_78354в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
└2╜
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_78359
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_78364в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╥2╧
(__inference_conv2d_3_layer_call_fn_78373в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
э2ъ
C__inference_conv2d_3_layer_call_and_return_conditional_losses_78384в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╥2╧
(__inference_conv2d_4_layer_call_fn_78393в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
э2ъ
C__inference_conv2d_4_layer_call_and_return_conditional_losses_78404в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
К2З
/__inference_max_pooling2d_3_layer_call_fn_78409
/__inference_max_pooling2d_3_layer_call_fn_78414в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
└2╜
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_78419
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_78424в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╤2╬
'__inference_flatten_layer_call_fn_78429в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ь2щ
B__inference_flatten_layer_call_and_return_conditional_losses_78435в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╧2╠
%__inference_dense_layer_call_fn_78450в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ъ2ч
@__inference_dense_layer_call_and_return_conditional_losses_78467в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╤2╬
'__inference_dense_1_layer_call_fn_78482в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ь2щ
B__inference_dense_1_layer_call_and_return_conditional_losses_78499в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
М2Й
'__inference_dropout_layer_call_fn_78504
'__inference_dropout_layer_call_fn_78509┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
┬2┐
B__inference_dropout_layer_call_and_return_conditional_losses_78514
B__inference_dropout_layer_call_and_return_conditional_losses_78526┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╤2╬
'__inference_dense_2_layer_call_fn_78541в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ь2щ
B__inference_dense_2_layer_call_and_return_conditional_losses_78558в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
╤2╬
'__inference_dense_3_layer_call_fn_78573в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ь2щ
B__inference_dense_3_layer_call_and_return_conditional_losses_78590в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
Р2Н
)__inference_dropout_1_layer_call_fn_78595
)__inference_dropout_1_layer_call_fn_78600┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╞2├
D__inference_dropout_1_layer_call_and_return_conditional_losses_78605
D__inference_dropout_1_layer_call_and_return_conditional_losses_78617┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p 

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╨2═
&__inference_output_layer_call_fn_78626в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
ы2ш
A__inference_output_layer_call_and_return_conditional_losses_78637в
Щ▓Х
FullArgSpec
argsЪ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
▓2п
__inference_loss_fn_0_78648П
З▓Г
FullArgSpec
argsЪ 
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *в 
▓2п
__inference_loss_fn_1_78659П
З▓Г
FullArgSpec
argsЪ 
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *в 
▓2п
__inference_loss_fn_2_78670П
З▓Г
FullArgSpec
argsЪ 
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *в 
▓2п
__inference_loss_fn_3_78681П
З▓Г
FullArgSpec
argsЪ 
varargs
 
varkw
 
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *в 
╙B╨
#__inference_signature_wrapper_77313sequential_input"Ф
Н▓Й
FullArgSpec
argsЪ 
varargs
 
varkwjkwargs
defaults
 

kwonlyargsЪ 
kwonlydefaults
 
annotationsк *
 
Ф2С
+__inference_random_flip_layer_call_fn_78686
+__inference_random_flip_layer_call_fn_78693┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╩2╟
F__inference_random_flip_layer_call_and_return_conditional_losses_78697
F__inference_random_flip_layer_call_and_return_conditional_losses_78755┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
Ь2Щ
/__inference_random_rotation_layer_call_fn_78760
/__inference_random_rotation_layer_call_fn_78767┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╥2╧
J__inference_random_rotation_layer_call_and_return_conditional_losses_78771
J__inference_random_rotation_layer_call_and_return_conditional_losses_78893┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
Ф2С
+__inference_random_zoom_layer_call_fn_78898
+__inference_random_zoom_layer_call_fn_78905┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 
╩2╟
F__inference_random_zoom_layer_call_and_return_conditional_losses_78909
F__inference_random_zoom_layer_call_and_return_conditional_losses_79023┤
л▓з
FullArgSpec)
args!Ъ
jself
jinputs

jtraining
varargs
 
varkw
 
defaultsЪ
p

kwonlyargsЪ 
kwonlydefaultsк 
annotationsк *
 п
 __inference__wrapped_model_75802К%&/09:CDIJWX]^ghmnwxAв>
7в4
2К/
sequential_input         &
к "/к,
*
output К
output         │
C__inference_conv2d_1_layer_call_and_return_conditional_losses_78304l/07в4
-в*
(К%
inputs         
к "-в*
#К 
0         
Ъ Л
(__inference_conv2d_1_layer_call_fn_78293_/07в4
-в*
(К%
inputs         
к " К         │
C__inference_conv2d_2_layer_call_and_return_conditional_losses_78344l9:7в4
-в*
(К%
inputs         	
к "-в*
#К 
0         	 
Ъ Л
(__inference_conv2d_2_layer_call_fn_78333_9:7в4
-в*
(К%
inputs         	
к " К         	 │
C__inference_conv2d_3_layer_call_and_return_conditional_losses_78384lCD7в4
-в*
(К%
inputs          
к "-в*
#К 
0         @
Ъ Л
(__inference_conv2d_3_layer_call_fn_78373_CD7в4
-в*
(К%
inputs          
к " К         @│
C__inference_conv2d_4_layer_call_and_return_conditional_losses_78404lIJ7в4
-в*
(К%
inputs         @
к "-в*
#К 
0         P
Ъ Л
(__inference_conv2d_4_layer_call_fn_78393_IJ7в4
-в*
(К%
inputs         @
к " К         P▒
A__inference_conv2d_layer_call_and_return_conditional_losses_78264l%&7в4
-в*
(К%
inputs         &
к "-в*
#К 
0         &
Ъ Й
&__inference_conv2d_layer_call_fn_78253_%&7в4
-в*
(К%
inputs         &
к " К         &д
B__inference_dense_1_layer_call_and_return_conditional_losses_78499^]^0в-
&в#
!К
inputs         ю
к "&в#
К
0         м
Ъ |
'__inference_dense_1_layer_call_fn_78482Q]^0в-
&в#
!К
inputs         ю
к "К         мг
B__inference_dense_2_layer_call_and_return_conditional_losses_78558]gh0в-
&в#
!К
inputs         м
к "%в"
К
0         d
Ъ {
'__inference_dense_2_layer_call_fn_78541Pgh0в-
&в#
!К
inputs         м
к "К         dв
B__inference_dense_3_layer_call_and_return_conditional_losses_78590\mn/в,
%в"
 К
inputs         d
к "%в"
К
0         2
Ъ z
'__inference_dense_3_layer_call_fn_78573Omn/в,
%в"
 К
inputs         d
к "К         2в
@__inference_dense_layer_call_and_return_conditional_losses_78467^WX0в-
&в#
!К
inputs         а
к "&в#
К
0         ю
Ъ z
%__inference_dense_layer_call_fn_78450QWX0в-
&в#
!К
inputs         а
к "К         юд
D__inference_dropout_1_layer_call_and_return_conditional_losses_78605\3в0
)в&
 К
inputs         2
p 
к "%в"
К
0         2
Ъ д
D__inference_dropout_1_layer_call_and_return_conditional_losses_78617\3в0
)в&
 К
inputs         2
p
к "%в"
К
0         2
Ъ |
)__inference_dropout_1_layer_call_fn_78595O3в0
)в&
 К
inputs         2
p 
к "К         2|
)__inference_dropout_1_layer_call_fn_78600O3в0
)в&
 К
inputs         2
p
к "К         2д
B__inference_dropout_layer_call_and_return_conditional_losses_78514^4в1
*в'
!К
inputs         м
p 
к "&в#
К
0         м
Ъ д
B__inference_dropout_layer_call_and_return_conditional_losses_78526^4в1
*в'
!К
inputs         м
p
к "&в#
К
0         м
Ъ |
'__inference_dropout_layer_call_fn_78504Q4в1
*в'
!К
inputs         м
p 
к "К         м|
'__inference_dropout_layer_call_fn_78509Q4в1
*в'
!К
inputs         м
p
к "К         мз
B__inference_flatten_layer_call_and_return_conditional_losses_78435a7в4
-в*
(К%
inputs         P
к "&в#
К
0         а
Ъ 
'__inference_flatten_layer_call_fn_78429T7в4
-в*
(К%
inputs         P
к "К         а:
__inference_loss_fn_0_78648Wв

в 
к "К :
__inference_loss_fn_1_78659]в

в 
к "К :
__inference_loss_fn_2_78670gв

в 
к "К :
__inference_loss_fn_3_78681mв

в 
к "К э
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_78319ЮRвO
HвE
CК@
inputs4                                    
к "HвE
>К;
04                                    
Ъ ╢
J__inference_max_pooling2d_1_layer_call_and_return_conditional_losses_78324h7в4
-в*
(К%
inputs         
к "-в*
#К 
0         	
Ъ ┼
/__inference_max_pooling2d_1_layer_call_fn_78309СRвO
HвE
CК@
inputs4                                    
к ";К84                                    О
/__inference_max_pooling2d_1_layer_call_fn_78314[7в4
-в*
(К%
inputs         
к " К         	э
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_78359ЮRвO
HвE
CК@
inputs4                                    
к "HвE
>К;
04                                    
Ъ ╢
J__inference_max_pooling2d_2_layer_call_and_return_conditional_losses_78364h7в4
-в*
(К%
inputs         	 
к "-в*
#К 
0          
Ъ ┼
/__inference_max_pooling2d_2_layer_call_fn_78349СRвO
HвE
CК@
inputs4                                    
к ";К84                                    О
/__inference_max_pooling2d_2_layer_call_fn_78354[7в4
-в*
(К%
inputs         	 
к " К          э
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_78419ЮRвO
HвE
CК@
inputs4                                    
к "HвE
>К;
04                                    
Ъ ╢
J__inference_max_pooling2d_3_layer_call_and_return_conditional_losses_78424h7в4
-в*
(К%
inputs         P
к "-в*
#К 
0         P
Ъ ┼
/__inference_max_pooling2d_3_layer_call_fn_78409СRвO
HвE
CК@
inputs4                                    
к ";К84                                    О
/__inference_max_pooling2d_3_layer_call_fn_78414[7в4
-в*
(К%
inputs         P
к " К         Pы
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_78279ЮRвO
HвE
CК@
inputs4                                    
к "HвE
>К;
04                                    
Ъ ┤
H__inference_max_pooling2d_layer_call_and_return_conditional_losses_78284h7в4
-в*
(К%
inputs         &
к "-в*
#К 
0         
Ъ ├
-__inference_max_pooling2d_layer_call_fn_78269СRвO
HвE
CК@
inputs4                                    
к ";К84                                    М
-__inference_max_pooling2d_layer_call_fn_78274[7в4
-в*
(К%
inputs         &
к " К         б
A__inference_output_layer_call_and_return_conditional_losses_78637\wx/в,
%в"
 К
inputs         2
к "%в"
К
0         
Ъ y
&__inference_output_layer_call_fn_78626Owx/в,
%в"
 К
inputs         2
к "К         ╢
F__inference_random_flip_layer_call_and_return_conditional_losses_78697l;в8
1в.
(К%
inputs         &
p 
к "-в*
#К 
0         &
Ъ ║
F__inference_random_flip_layer_call_and_return_conditional_losses_78755pў;в8
1в.
(К%
inputs         &
p
к "-в*
#К 
0         &
Ъ О
+__inference_random_flip_layer_call_fn_78686_;в8
1в.
(К%
inputs         &
p 
к " К         &Т
+__inference_random_flip_layer_call_fn_78693cў;в8
1в.
(К%
inputs         &
p
к " К         &║
J__inference_random_rotation_layer_call_and_return_conditional_losses_78771l;в8
1в.
(К%
inputs         &
p 
к "-в*
#К 
0         &
Ъ ╛
J__inference_random_rotation_layer_call_and_return_conditional_losses_78893p¤;в8
1в.
(К%
inputs         &
p
к "-в*
#К 
0         &
Ъ Т
/__inference_random_rotation_layer_call_fn_78760_;в8
1в.
(К%
inputs         &
p 
к " К         &Ц
/__inference_random_rotation_layer_call_fn_78767c¤;в8
1в.
(К%
inputs         &
p
к " К         &╢
F__inference_random_zoom_layer_call_and_return_conditional_losses_78909l;в8
1в.
(К%
inputs         &
p 
к "-в*
#К 
0         &
Ъ ║
F__inference_random_zoom_layer_call_and_return_conditional_losses_79023pГ;в8
1в.
(К%
inputs         &
p
к "-в*
#К 
0         &
Ъ О
+__inference_random_zoom_layer_call_fn_78898_;в8
1в.
(К%
inputs         &
p 
к " К         &Т
+__inference_random_zoom_layer_call_fn_78905cГ;в8
1в.
(К%
inputs         &
p
к " К         &▓
F__inference_rescaling_1_layer_call_and_return_conditional_losses_78244h7в4
-в*
(К%
inputs         &
к "-в*
#К 
0         &
Ъ К
+__inference_rescaling_1_layer_call_fn_78236[7в4
-в*
(К%
inputs         &
к " К         &╘
G__inference_sequential_1_layer_call_and_return_conditional_losses_77143И%&/09:CDIJWX]^ghmnwxIвF
?в<
2К/
sequential_input         &
p 

 
к "%в"
К
0         
Ъ ┌
G__inference_sequential_1_layer_call_and_return_conditional_losses_77236Оў¤Г%&/09:CDIJWX]^ghmnwxIвF
?в<
2К/
sequential_input         &
p

 
к "%в"
К
0         
Ъ ╔
G__inference_sequential_1_layer_call_and_return_conditional_losses_77519~%&/09:CDIJWX]^ghmnwx?в<
5в2
(К%
inputs         &
p 

 
к "%в"
К
0         
Ъ ╨
G__inference_sequential_1_layer_call_and_return_conditional_losses_77925Дў¤Г%&/09:CDIJWX]^ghmnwx?в<
5в2
(К%
inputs         &
p

 
к "%в"
К
0         
Ъ л
,__inference_sequential_1_layer_call_fn_76638{%&/09:CDIJWX]^ghmnwxIвF
?в<
2К/
sequential_input         &
p 

 
к "К         ▓
,__inference_sequential_1_layer_call_fn_77056Бў¤Г%&/09:CDIJWX]^ghmnwxIвF
?в<
2К/
sequential_input         &
p

 
к "К         б
,__inference_sequential_1_layer_call_fn_77358q%&/09:CDIJWX]^ghmnwx?в<
5в2
(К%
inputs         &
p 

 
к "К         з
,__inference_sequential_1_layer_call_fn_77409wў¤Г%&/09:CDIJWX]^ghmnwx?в<
5в2
(К%
inputs         &
p

 
к "К         ─
E__inference_sequential_layer_call_and_return_conditional_losses_76211{JвG
@в=
3К0
random_flip_input         &
p 

 
к "-в*
#К 
0         &
Ъ ═
E__inference_sequential_layer_call_and_return_conditional_losses_76224Гў¤ГJвG
@в=
3К0
random_flip_input         &
p

 
к "-в*
#К 
0         &
Ъ ╣
E__inference_sequential_layer_call_and_return_conditional_losses_77945p?в<
5в2
(К%
inputs         &
p 

 
к "-в*
#К 
0         &
Ъ ┴
E__inference_sequential_layer_call_and_return_conditional_losses_78231xў¤Г?в<
5в2
(К%
inputs         &
p

 
к "-в*
#К 
0         &
Ъ Ь
*__inference_sequential_layer_call_fn_75831nJвG
@в=
3К0
random_flip_input         &
p 

 
к " К         &д
*__inference_sequential_layer_call_fn_76204vў¤ГJвG
@в=
3К0
random_flip_input         &
p

 
к " К         &С
*__inference_sequential_layer_call_fn_77930c?в<
5в2
(К%
inputs         &
p 

 
к " К         &Щ
*__inference_sequential_layer_call_fn_77941kў¤Г?в<
5в2
(К%
inputs         &
p

 
к " К         &╞
#__inference_signature_wrapper_77313Ю%&/09:CDIJWX]^ghmnwxUвR
в 
KкH
F
sequential_input2К/
sequential_input         &"/к,
*
output К
output         