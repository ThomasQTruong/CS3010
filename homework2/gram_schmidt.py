"""
  Implementation of the Gram-Schmidt process.
"""

import math

def magnitude(vector):
  """
    Calculates the magnitude of the vector.
    
    Args:
      vector - the vector to calculate the magnitude for.
    
    Returns:
      The magnitude of the vector.
  """
  total = 0
  for num in vector:
    total += num ** 2

  return math.sqrt(total)


def subtract_vectors(vector1, vector2):
  """
    Subtracts vector1 and vector2.
    
    Args:
      vector1 - the first vector.
      vector2 - the second vector.
    
    Returns:
      The resulting vector after subtraction.
  """
  subtracted = []
  for i in range(len(vector1)):
    subtracted.append(vector1[i] - vector2[i])

  return subtracted


def multiply_row_col_vectors(row_vector, col_vector):
  """
    Multiplies 2 vectors together, where one is transposed (a^T * b).
    
    Args:
      row_vector - vector in row format.
      col_vector - vector in column format.
    
    Returns:
      The sum of the two vectors multiplied (summation(a_i^T * b_i)).
  """
  total = 0
  for i in range(len(row_vector)):
    total += row_vector[i] * col_vector[i]

  return total


def scale_vector(vector, constant):
  """
    Scales vector by a constant.

    Args:
      vector - the vector to get multiplied.
      constant - the number to multiply by.
    
    Returns:
      A vector that was scaled by the constant.
  """
  scaled_vector = []
  for i in range(len(vector)):
    scaled_vector.append(vector[i] * constant)

  return scaled_vector


def projection(v, u):
  """
    Implementing the orthogonal projection formula:
      w = proj_v(u) = [(u^Tv) / (v^Tv)]v
        
    Args:
      v - vector to project onto u.
      u - vector to get projected on.
    Returns:
      The orthogonal projection of v onto u.
  """
  # v * [(u^Tv) / (v^Tv)]
  w = scale_vector(v,
        float(multiply_row_col_vectors(u, v)) / multiply_row_col_vectors(v, v))
  return w


def gram_schmidt(v_list):
  """
      Orthonormalizes a set of vectors.
      
      Args:
          v_list - the set of vectors.
      
      Returns:
          The set of orthonormalizes vectors.
  """
  e = []  # Orthonormalized vectors.
  u = []
  # For every vector in the list.
  for i in range(len(v_list)):
    # Add current vector into u.
    u.append(v_list[i])
    # For every existing u.
    for j in range(i):
      u[i] = subtract_vectors(u[i], projection(u[j], v_list[i]))
    e.append(scale_vector(u[i], 1/magnitude(u[i])))

  return e


if __name__ == "__main__":
  vector_list = [[2, 1, -1],
                 [1, 1, 0],
                 [0, -1, 1]]
  results = gram_schmidt(vector_list)

  print("Orthonormalized Vectors:")
  for result in results:
    print(result)
