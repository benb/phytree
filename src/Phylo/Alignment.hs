module Phylo.Alignment where
import Phylo.Alignment.Parsers
import Phylo.Tree hiding (names)
import qualified Phylo.Tree as Tree
import qualified Data.ByteString.Lazy.Char8 as L
import Control.Monad
import Data.List
import Data.Char
import Debug.Trace
import qualified Data.HashMap as HM
import Text.JSON
import Data.Maybe



parseAlignmentString :: Monad m =>  ([L.ByteString] -> m [(String,String)]) -> L.ByteString -> m ListAlignment

parseAlignmentString parser input = (liftM2 safeListAlignment names seqs) where 
                                mydata = (liftM $ sortBy sortX) (parser (L.lines input))
                                sortX (a,b) (c,d) = compare a c
                                names = (liftM $ map fst) mydata
                                seqs = (liftM $ map snd) mydata 


parseAlignmentFile :: Monad m => ([L.ByteString] -> m [(String,String)]) -> String -> IO (m ListAlignment)
parseAlignmentFile parser name = parseAlignmentString parser `liftM` (L.readFile name)

parseUniversal :: Monad m => [L.ByteString] -> m [(String,String)]
parseUniversal [] = error "Empty alignment file"
parseUniversal (x:xs) = case x of 
                              str | str==L.empty -> parseUniversal xs
                                  | L.isPrefixOf (L.pack ">") str -> parseFasta (x:xs)
                                  | otherwise -> parsePhylip(x:xs)


                                

dropGaps :: ListAlignment -> [(String,String)]
dropGaps a = zip (names a) (map dropGap $ sequences a)

dropGap :: String -> String
dropGap xs = filter (not . isGapChar) xs

isGapChar :: Char -> Bool
isGapChar '-' = True
isGapChar '.' = True
isGapChar '~' = True
isGapChar x = False


type Name = String
type Sequence = [Char]
type Column = [Char]

bsort :: Ord a => [a] -> [a]
bsort s = case _bsort s of
               t | t == s    -> t
                 | otherwise -> bsort t
          where _bsort (x:x2:xs) | x > x2    = x2:(_bsort (x:xs))
                                 | otherwise = x:(_bsort (x2:xs))
                _bsort s = s

data NumberedColumn = NumberedColumn {coldata::[(Char,Int)]} deriving (Eq)

--This is not a real Ord
--It will only work with neighbour-neighbour comparisons from the starting alignment
--Hence - use Bubble Sort
--A real Ord impl would require each pair of columns to be examined in the 
--Context of the whole alignment

instance Ord NumberedColumn where 
                       compare (NumberedColumn x) (NumberedColumn y) = compare' x y Nothing where
                                --maintain sequence ordering
                                compare' ((x,i):xs) ((y,j):ys) ans | (not $ isGapChar x) && (not $ isGapChar y) = compare i j

                                -- two all gap columns!?
                                compare' [] [] Nothing = EQ 

                                -- all pairs are gap-base or base-gap --> arbitrary
                                compare' [] [] (Just ans) = ans

                                --first gap on left side -> set arbitrary answer and keep looking for base-base pairs
                                compare' ((gap,i):xs) ((y,j):ys) Nothing | (isGapChar gap) && (not $ isGapChar y) = compare' xs ys (Just GT)
                                compare' ((x,i):xs) ((gap,j):ys) Nothing | (isGapChar gap) && (not $ isGapChar x) = compare' xs ys (Just LT)

                                --existing ordering
                                compare' (x:xs) (y:ys) ans = compare' xs ys ans

sortAlignment (ListAlignment names seqs cols) = ListAlignment names (transpose ans) ans where
                                                  numbers = transpose $ numberifyBasic $ ListAlignment names seqs cols
                                                  numbCols = map NumberedColumn $ map (\(a,b)-> zip a b) $ zip cols numbers
                                                  reordered = bsort numbCols
                                                  ans = map (map fst) (map coldata reordered)

sortAlignmentBySeq (ListAlignment names seqs cols) = sortByName $ sortAlignment aln2 where
                                                                aln2 = ListAlignment newnames newseqs newcols 
                                                                (newnames,newseqs) = unzip $ sortBy (\(a,b) (c,d) -> compare (removeGaps b) (removeGaps d)) (zip names seqs)
                                                                removeGaps (x:xs) | isGapChar x = removeGaps xs
                                                                                  | otherwise = x:(removeGaps xs)
                                                                removeGaps [] = []
                                                                newcols = transpose newseqs

sortByName (ListAlignment names seqs cols) = ListAlignment newnames newseqs (transpose newseqs) where
                                                 (newnames,newseqs) = unzip $ sortBy (\(a,b) (c,d) -> compare a c) $ zip names seqs 

                                                        



gapPos :: Sequence -> [(Int,Int)]
gapPos s = gapPos' s [] Nothing 0

absoluteGapPos :: [(Int,Int)] -> [(Int,Int)]
absoluteGapPos s = map firstTwo (absoluteGapPos' s) where
                        firstTwo (a,b,c) = (a,b)

absoluteGapPos' :: [(Int,Int)] -> [(Int,Int,Int)]
absoluteGapPos' [] = [] 
absoluteGapPos' ((i,j):[]) = (i,j,0):[]
absoluteGapPos' ((i,j):xs) = (i+myoffset,j,myoffset) : agpTail where
                                agpTail = absoluteGapPos' xs
                                offset ((a,b,c):ys) = c + b
                                myoffset = offset agpTail


gapPos' :: Sequence -> [(Int,Int)] -> (Maybe Int) -> Int -> [(Int,Int)]
--end of sequence
gapPos' [] list Nothing pos = list 
gapPos' [] list (Just i) pos = (pos,i):list
--open a gap
gapPos' (gap:xs) list Nothing pos | isGapChar gap = gapPos' xs list (Just 1) (pos) 
--extend a gap
gapPos' (gap:xs) list (Just i) pos | isGapChar gap = gapPos' xs list (Just (i+1)) (pos)
--close a gap
gapPos' (x:xs) list (Just i) pos = gapPos' xs ((pos,i):list) Nothing (pos+1)
--no gap
gapPos' (x:xs) list Nothing pos = gapPos' xs list Nothing (pos+1)



data ListAlignment = ListAlignment {names ::  [Name],
                            sequences :: [Sequence],
                            columns :: [Column]} deriving Show

instance Eq ListAlignment where 
        (==) a b = (names a) == (names b) && (sequences a) == (sequences b) --assume cols are ok

instance JSON ListAlignment where 
        showJSON (ListAlignment names seqs cols) = showJSON ("Alignment",names,seqs)
        readJSON (JSArray [(JSString tag),(JSArray seqs),(JSArray vals)]) = Error "unimplemented" -- quickListAlignment names seqs


quickListAlignment :: [Name] -> [Sequence] -> ListAlignment
quickListAlignment names sequences = ListAlignment names sequences (transpose sequences)

saneAlignment (ListAlignment names sequences columns) = saneAlignment' names sequences columns
saneAlignment' names sequences columns = hasDuplicate names
hasDuplicate (x:y:z) | x==y = Just $ "Duplicate sequence identifier " ++ x
                     | otherwise = hasDuplicate (y:z)
hasDuplicate z = Nothing

safeListAlignment names sequences = case (saneAlignment ans) of 
                                        Just errmsg -> error errmsg
                                        Nothing -> ans
                                        where ans = removeAllGaps $ quickListAlignment names (map (map toUpper) sequences)

fromColumnListAlignment :: [Name] -> [Column] -> ListAlignment
fromColumnListAlignment names cols = ListAlignment names (transpose cols) cols


toFasta :: ListAlignment -> [String]
toFasta aln = stringList where --foldl (++) "" stringList where 
                 stringList = map toSeqStr seqList
                 seqList = zip (names aln) (sequences aln)
                 toSeqStr :: (String,String) -> String
                 toSeqStr (name,seq) = ">" ++ name ++ "\n" ++ seq ++ "\n"

removeAllGaps :: ListAlignment -> ListAlignment
removeAllGaps (ListAlignment names seqs cols) = fromColumnListAlignment names (removeAllGaps' cols)

removeAllGaps' :: [Column] -> [Column]
removeAllGaps' = filter notAllGap where 
 
notAllGap :: Column -> Bool
notAllGap (gap:[]) | isGapChar gap = False
notAllGap (gap:xs) | isGapChar gap = notAllGap xs
notAllGap (x:xs) = True

hasGap :: Column -> Bool
hasGap [] = False
hasGap (gap:xs) | isGapChar gap = True
hasGap (x:xs) = hasGap xs


--orderGaps :: ListAlignment -> ListAlignment
--orderGaps (ListAlignment names seqs cols) = fromColumnListAlignment names (orderGaps' cols)


--orderGaps' :: [Column] -> [Column] -> [Column]
--orderGaps' [] x = x:[]
--orderGaps' (x:xs) [] | hasGap x = orderGaps' xs x:[]
--orderGaps' (x:xs) [] = orderGaps' xs []
--orderGaps (x:xs) (y:ys) | canPushTogether x y 
                        
numberifyBasic :: ListAlignment -> [[Int]]
numberifyBasic aln = map nfy myseqs where
        myseqs = sequences aln
        nfy = numberMap 0
        numberMap i [] = []
        numberMap i (gap:xs) |isGapChar gap = -1 : numberMap i xs
        numberMap i (x:xs) = i : numberMap (i+1) xs

numberifyGap :: ListAlignment -> [[Int]]
numberifyGap aln = map nfy myseqs where
        myseqs = sequences aln
        nfy = numberMap 0
        numberMap i [] = []
        numberMap i (gap:xs) | isGapChar gap = (-(i+1)) : numberMap i xs
        numberMap i (x:xs) = i : numberMap (i+1) xs

numberifyGapTree :: Node -> ListAlignment -> [[(Int,Maybe Node)]]
numberifyGapTree tree aln = transpose $ nfy (columns aln) where
        nfy :: [Column]  -> [[(Int, Maybe Node)]]
        nfy colList = numberMap (map (\x->0) (head (columns aln))) colList
        numberMap :: [Int] -> [Column] -> [[(Int, Maybe Node)]]
        numberMap y [] = []
        numberMap y (x:xs) = (snd ans) : (numberMap (fst ans) xs) where
                              ans = numberMap' y x $ names aln
                              gapNums = unrootedSplitsFor tree gapNames
                              gapNames = map (\x-> fst x) $ filter (\t -> isGapChar (snd t)) $ zip (names aln) x
                              numberMap':: [Int] -> Column -> [String] -> ([Int],[(Int,Maybe Node)]) 
                              numberMap' [] [] [] = ([],[])

                              numberMap' (a:as) (gap:bs) (name:cs) | isGapChar gap = (a:(fst ans2),((-a-1,Just $ getNode name)):(snd ans2)) where
                                                                                             ans2 = numberMap' as bs cs
                              numberMap' (a:as) (b:bs) (name:cs) = (a+1:(fst ans2),(a,Nothing):(snd ans2)) where
                                                                                             ans2 = numberMap' as bs cs
                              getNode::String -> Node
                              getNode name = case (HM.lookup name gapNums) of 
                                Nothing -> error $ "Can't find gap for " ++ name
                                Just a -> a

compatible :: Node -> ListAlignment -> Bool
compatible tree aln = (sort $ Tree.names tree) == names aln

-- Incompat describes an incompatability between two alignments
-- The boolean value describes whether an incompatability 
-- renders comparison impossible (i.e. different symbols can be ignored to make
-- a comparison, different lengths cannot ) 
newtype Incompat = Incompat (Bool,String) 

incompatibilities :: ListAlignment -> ListAlignment -> [Incompat]
incompatibilities  (ListAlignment namesx seqsx colsx) (ListAlignment namesy seqsy colsy) | (length namesx) > (length namesy) = (Incompat (True,"incompatible alignments: first alignment has " ++ (show $ length $ namesx) ++ " sequences, second has " ++ (show $ length $ namesy))) : (incompatibilities' namesx seqsx namesy seqsy [])
incompatibilities (ListAlignment namesx seqsx colsx) (ListAlignment namesy seqsy colsy) = incompatibilities' namesx seqsx namesy seqsy [] 
incompatibilities' (namex:namexs) (seqx:seqxs) (namey:nameys) (seqy:seqys) incompat = incompatibilities' namexs seqxs nameys seqys $ compatibleSeqs namex seqx namey seqy incompat

incompatibilities' [] [] [] [] incompat = incompat 

compatibleSeqs namex seqx namey seqy incompat | namex /= namey = (Incompat(True,"incompatible alignments: " ++ " incompatible sequence names " ++ namex ++" != " ++ namey)) : incompat
compatibleSeqs namex seqx namey seqy incompat = compatibleSeqs' namex (dropGap seqx) (dropGap seqy) 0 incompat 

compatibleSeqs' :: String -> String -> String -> Int -> [Incompat] -> [Incompat]
compatibleSeqs' namex (x:seqx) (y:seqy) pos incompat | x/=y = compatibleSeqs' namex seqx seqy (pos+1) ((Incompat (False,"incompatible alignments: sequence " ++ namex ++ " differs between two alignments " ++ difference)):incompat)  where
                                                                          difference = " as there is a different character at ungapped position " ++ (show pos) ++ " ( " ++ (show x) ++ " vs " ++ (show y) ++ ")" 
compatibleSeqs' namex (x:seqx) (y:seqy) pos incompat = compatibleSeqs' namex seqx seqy (pos+1) incompat
compatibleSeqs' namex (x:seqx) [] pos incompat = (Incompat (True,"incompatible alignments: sequence " ++ namex ++ " differs between two alignments as they have different lengths ( " ++ (show (pos + (length seqx) + 1)) ++ " vs " ++ (show pos) ++ ")"))  : incompat
compatibleSeqs' namex [] (y:seqy) pos incompat = (Incompat (True,"incompatible alignments: sequence " ++ namex ++ " differs between two alignments as they have different lengths ( " ++ (show pos) ++ " vs " ++ (show (pos + (length seqy) + 1)) ++ ")"))  : incompat
compatibleSeqs' namex [] [] pos incompat = incompat


scaledAAFrequencies = scaledFrequencies "ARNDCQEGHILKMFPSTWYV"
scaledFrequencies :: [Char] -> ListAlignment -> [Double]
scaledFrequencies wanted aln = map (\i-> (fromIntegral i)/total) rawfreq where
                               total = fromIntegral $ foldr (+) 0 rawfreq
                               rawfreq = frequencies wanted aln


frequencies :: [Char] -> ListAlignment -> [Int]
frequencies wanted aln = ans where
                          hm = foldr (\k h -> HM.insert k 0 h) HM.empty wanted
                          hm2 = frequencies' hm (sequences aln)
                          ans = map (\k -> fromJust $ HM.lookup k hm2) wanted
frequencies' hm ((x:xs):xxs) = frequencies' hm2 (xs:xxs) where
                                hm2 = HM.update (\x -> Just $ x+1) x hm
frequencies' hm ([]:xxs) = frequencies' hm xxs
frequencies' hm ([]) = hm

