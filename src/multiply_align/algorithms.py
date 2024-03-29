import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from .nn_model import create_nn_score_dt
import os

DIR = os.path.dirname(__file__)


# ================================================================================
# Define an alignment between two primers
#
# ================================================================================


@dataclass(order=True)
class PrimerAlignment:
    """
    Represent the alignment of two primers
    """

    primer1: str = field(compare=False)
    primer2: str = field(compare=False)
    score: float = field(compare=True)
    alignment: str = field(compare=False, repr=False)
    # index: int=field(compare=False, repr=False)


# ================================================================================
# Abstract base class for various primer alignment algorithms
#
# ================================================================================


class AlignmentAlgorithm(ABC):
    """
    Alignment algorithm for a pair of primers

    """

    rc_map = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def __init__(self):
        pass

    def set_primers(self, primer1, primer2):
        """
        Set a pair for primers to align

        """
        self.primer1 = primer1
        self.primer2 = primer2

        self.score = None  # reset

    @abstractmethod
    def load_parameters():
        pass

    @abstractmethod
    def align():
        """
        Align the primers

        """
        pass

    def get_primer_alignment(self):
        """
        Return an alignment object

        """
        return PrimerAlignment(
            primer1=self.primer1,
            primer2=self.primer2,
            score=self.score,
        )


# ================================================================================
# Concrete primer alignment algorithms
#
# ================================================================================


class PrimerDimerLike(AlignmentAlgorithm):
    """
    Align two primers using an algorithm like the one described
    by Johnston et al. (2019) Sci Reports
    
    Primary idea is to only allow:
    - Ungapped alignments
    - With 5' overhangs (i.e. extensible)
    
    And then to add a bonus if either 3' end is complementary in
    the highest scoring alignment
    
    
    """
    
    param_path = os.path.join(DIR, "parameters.json")
    
    def load_parameters(self):
        """
        Load parameters necessary for Primer Dimer algorithm,
        and set as attributes
        
        """
        # Load parameter JSON
        params = {    
                    "end_length": 4,
                    "end_bonus": -0.5,
                    "double_mismatch_score": 0.2
                }

        # Load nearest neighbour model
        self.nn_scores = create_nn_score_dt(
            match_json=os.path.join(DIR, "match.json"),  # this is a path
            single_mismatch_json=os.path.join(DIR, "single_mismatch.json"),  # this is a path
            double_mismatch_score=params['double_mismatch_score']  # this is a float
        )

        # Load penalties
        self.end_length = params["end_length"]
        #self.end_penalty = params["end_penalty"]
        self.end_bonus = params["end_bonus"]
        
    @staticmethod
    def _calc_linear_extension_bonus(matching, 
                                     overhang_left, 
                                     overhang_right, 
                                     end_length, 
                                     end_bonus):
        """
        Calculate a bonus score for primer-dimer alignments that would allow for *extension*

        params
            matching : list[bool]
                List of booleans indicating whether or not bases matched
                for this alignment position
            overhang_left : bool
                Is there an overhang on the left side of the matches;
                i.e. would extension be possible?
            overhang_right : bool
                As above, but right side.
            end_length : int
                Number of bases to consider for end bonus.
            end_bonus : float
                Bonus to add per aligned, extendible base.

        returns
            _ : float
                Bonus score in [-2*end_length*end_bonus, 0].

        """

        # Left end
        left_end = 0
        for match in matching[:end_length]:
            if not overhang_left:
                break
            if not match:
                break
            left_end += 1

        # Right end
        right_end = 0
        for match in matching[::-1][:end_length]:
            if not overhang_right:
                break
            if not match:
                break
            right_end += 1

        return (left_end + right_end) * end_bonus
        
    def align(self):
        """
        Align primers;
        
        Finding highest score and its associated
        start position
        
        """
        
        # Identify longer and shorter primer
        primers = [self.primer1, self.primer2]
        primers.sort(key=len, reverse=True)
        l, s = primers
        s = s[::-1]
        nL, nS = len(l), len(s)

        # Iterate over each start position
        best_start = 0
        best_score = 10
        best_matching = []
        for i in range(nL - 1):

            # Compute score for alignment
            matching = []
            current_score = 0
            penalties = 0
            for j in range(nS - 1):

                # Compute psuedo-Gibb's
                l_bases = l[(i+j):(i+j+2)]
                s_bases = s[j:(j+2)]
                nn = f"{l_bases}/{s_bases}"
                current_score += self.nn_scores[nn]

                # Compute if matching, for left-end
                matching.append(self.rc_map[s[j]] == l[i+j])

                # Stop if reached last dinucleotide of longer primer
                if i + j == nL - 2:
                    break

            # Add matching state for last nucleotide
            matching.append(self.rc_map[s[j+1]] == l[i+j+1])

            # Compute end bonus
            current_score += self._calc_linear_extension_bonus(
                matching,
                overhang_left=i > 0,
                overhang_right=len(matching) < nS, 
                end_length=self.end_length,
                end_bonus=self.end_bonus
            )

            # Update if this is the new best score
            if current_score <= best_score:
                best_score = current_score
                best_start = i
                best_matching = matching
                
        # Assign to instance variables
        self.score = best_score
        self.best_start = best_start
        
        # Helps with getting alignment string
        self.best_matching = best_matching
        self.s = s
        self.l = l
        
    