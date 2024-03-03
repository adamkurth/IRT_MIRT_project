## n - Parameter Logistic Models (nPL)

For a through understanding of MIRT and Item Response Theory (IRT), I want to describe each logistic model. I generalized to *nPL* for \( n = \{ 1,2,3 \}\) parameters.

## Parameters:

#### Ability \( (\theta) \):

- In IRT, the *ability* parameter is denoted as \(\theta\). This refers to the latent trait or *characteristic being measured by the test*. This could be the acedemic proficiciency in a subject to a personality trait. 
  
- (e.g.) In a math test, \( \theta \) represents the student's mathematical ability. A higher \( \theta \) indicates a higher level in math.

#### Difficulty \( (b)\):

- The *difficulty* parameter \( b \), indicates how challenging the item is. Each item in the test has this difficulty level.

- For 1PL, this is crucial this parameter is crucial. A higher \( b \) value means that the question is more difficult. 
- (e.g.) In the math test, a question on basic addition would have a low \( b \) value (easy), while a more complex calculus question has a higher \( b \) value. 

#### Discrimination \( (a) \): 

- This *discrimination* parameter is introduced in the 2PL, denoted \( a \), indicates how well the item can differentiate between the individuals of diffferent abilities. 

- A higher \( a \) value means that the item is better at distinguishing between those who have mastered the material and those who haven't.

- (e.g.) A well designed quesition that most students with a good understanding of algebra get right, but those who get it wrong, would have a high value. A quesiton that has both high and low ability tend to answer similarity (either correct/incorrect) would have a low value.

### Guessing \( (c) \):

- The *guessing* parameter is unique to the 3PL model. It represents the probability that a person with low ability could guess the correct answer.

- \( c \) is particularily relavent in multiple-choice tests where random guesses could lead to correct answers. 

- (e.g.) In the multiple-choice math test, with four options per question, a question where one could easily eliminate two wrong answers can have a higher \( c \) value since the probability of guessing between two options is higher. 


### One Parameter Logistic Model: *Rasch Model*

- Equation: \[ P(\theta) = \frac{1}{1+e^{-a(\theta - b)}} \]

- \( \theta \): Ability of the individual. 

- \( b \): Difficulty of the item.

- The 1PL, the only item that is being estimated is \( b \). It's assumed that all items have the same discrimination power (a = 1), and there's no guessing (c = 0).

- This model posits the probability of a correct response is a logistic function of the difference between the person's ability and the item's difficulty.

### Two Parameter Logistic Model: *(2PL)*

- Adds the *discrimination* parameter.

- Equation: \[ P(\theta) = c + (1-c)\frac{1}{1+e^{-a(\theta - b)}} \]

- \( \theta \): Ability of the individual. 

- \( a \): Discrimination of the item.

- \( b \): Difficulty of the item.

- Both the difficulty \( b \) and the discrimination \( a \) of an item are estimated. The discrimination parameter indicates how well an item distinguishes between the individuals with different levels of ability. 

- The higher the discrimination parameter means the item is better at differentiating between individuals with abilities above or below the item's difficulty level.

### Three Parameter Logistic Model: *(3PL)*

- Further extends the modelk by including a guessing parameter. 

- Equation: \[ P(\theta) = c + (1-c)\frac{1}{1+e^{-a(\theta - b)}} \]

- \( \theta \): Ability of the individual. 

- \( a \): Discrimination of the item.

- \( b \): Difficulty of the item.

- \( c \): Guessing parameter. 

- In the 3PL model, it estimates the difficulty \( b \), discrimination \( a \), and guessing parameters \( c \). \( c \) accounts for the likelihood that a low-ability \( \theta \) individual might guess the answer correctly.
- Model is often used in multiple choice tests where the possibility of guessing correctly is non-zero.
  
\[ P(\theta) = c + (1-c)\frac{1}{1+e^{-a(\theta - b)}} \]
