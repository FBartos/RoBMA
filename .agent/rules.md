# RoBMA R Code Style Rules

When writing or modifying R code in this package, strictly follow these rules:

## 1. Vertical Alignment
You MUST maintain vertical alignment for readability in the following code blocks:

### Assignments
Align the assignment arrow `<-` when multiple assignments occur in a sequence.
```r
# Correct
is_multilevel     <- .is_multilevel(x)
is_weightfunction <- .is_weightfunction(x)

# Incorrect
is_multilevel <- .is_multilevel(x)
is_weightfunction <- .is_weightfunction(x)
```

### Function Arguments
Align arguments when a function call spans multiple lines.
```r
# Correct
result <- my_function(
  first_arg  = value1,
  second_arg = value2
)
```

### Return Lists
Align the elements when returning a list.
```r
# Correct
return(list(
  val1 = x,
  val2 = y
))
```

## 2. Spacing and Indentation
- Use **2 spaces** for indentation.
- Always place spaces around assignment operators `<-` and other binary operators.
- Always place a space after a comma.
- Leave an empty line immediately after the opening brace `{` of a function definition.
