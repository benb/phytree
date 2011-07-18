module Phylo.Alignment.Parsers where
import qualified Data.ByteString.Lazy.Char8 as L
import Data.List
import Data.Char
import Control.Applicative
import Debug.Trace

--This parser has been briefly tested to work with the "relaxed" phylip format
--It is intended to handle interleaved and sequential versions, but relies on spaces between 
--The sequence name and the sequence data
--It also won't work if sequence names are repeated
parsePhylip :: Monad m =>  [L.ByteString] -> m [(String,String)]
parsePhylip = parsePhylipHeader
parsePhylipHeader :: Monad m => [L.ByteString] -> m [(String,String)]


parsePhylipHeader (x:xs)  = (case (ans,nChar) of
                                (Just a,Just c) | find (\x-> length (snd x)/=c) a == Nothing -> return a
                                _ -> fail "Can't parse alignment") where
                                header = L.dropWhile (==' ') x
                                t1 = L.readInt header
                                nTaxa = fst <$> t1
                                t2 = snd <$> t1 >>= L.readInt . L.dropWhile isSpace
                                nChar = fst <$> t2
                                filteredxs = filter (\x->L.empty /= (L.filter (not . isSpace) x)) xs
                                ans = ((parseBeginPhylipBody filteredxs) <$> nTaxa <*> nChar)


deByteString :: [(L.ByteString,L.ByteString)] -> [(String,String)]
deByteString = map (\(x,y) -> (L.unpack x,L.unpack$ L.filter (not . isSpace ) y))
parseBeginPhylipBody :: [L.ByteString] -> Int -> Int -> [(String,String)]
{-parseBeginPhylipBody (line:lines) nTaxa nChar | trace "OK1" False  = map (\(x,y) -> (L.unpack x,L.unpack y)) (parsePamlBody (line:lines) nTaxa nChar)-}
parseBeginPhylipBody (line:lines) nTaxa nChar | isPaml line = deByteString $ parsePamlBody (line:lines) nTaxa nChar
parseBeginPhylipBody (line:line2:lines) nTaxa nChar | isSequential (line,line2) = deByteString $ parseSeqPhylipBody (line:line2:lines) nTaxa nChar 
parseBeginPhylipBody (line:line2:lines) nTaxa nChar | isInterleaved (line,line2) = deByteString $ parseInterPhylipBody (line:line2:lines) nTaxa nChar 

isPaml line = 1 == (length $ (filter (not . L.null) (L.split ' ' line)))
isSequential (line1,line2) = isSpace $ L.head line2
isInterleaved (line1,line2) = not $ isSpace $ L.head line2

parsePamlBody :: [L.ByteString] -> Int -> Int -> [(L.ByteString,L.ByteString)]
{-parsePamlBody a b c | trace "OK" False = undefined-}
parsePamlBody lines 0 nChar  = [] 
parsePamlBody (line:lines) nTaxa nChar = sequence : parsePamlBody ltail (nTaxa-1) nChar where
                                       (name,remainder) = L.break isSpace line
                                       (seq,ltail) = getSeq nChar (remainder:lines) L.empty
                                       sequence = (name,seq)
parseSeqPhylipBody = parsePamlBody

getSeq :: Int -> [L.ByteString] -> L.ByteString -> (L.ByteString,[L.ByteString])
{-getSeq nChar lines seq | trace ((show nChar) ++ "\n" ++ (show (map L.unpack lines)) ++"\n" ++ (show $ L.unpack seq)) False = undefined-}
getSeq nChar lines seq | nChar < 0 = error "Failed to parse file"
getSeq 0 lines seq= (seq,lines)
getSeq nChar (line:xs) seq = getSeq (nChar- (fromIntegral len)) (xs) (seq `L.append` seq2) where
                                    seq2 = L.filter (not . isSpace) line
                                    len = L.length seq2
                                

{-parseInterPhylipBody lines nTaxa nChar | trace "Inter" False = undefined-}
parseInterPhylipBody lines nTaxa nChar = addLines lastLines $ firstInterPass firstLines where
                                                (firstLines,lastLines) = splitAt nTaxa lines
                                                firstInterPass =  map (L.break isSpace)
                                                addLines [] destination = destination
                                                addLines lines destination = addLines remainder $  map (\(seq2,(name,seq1)) -> (name,seq1 `L.append` seq2)) $ zip firstLines destination where
                                                                (firstLines,remainder) = splitAt nTaxa lines
                                                                                



--Fasta format
parseFasta' :: [(String,String)] -> [L.ByteString] -> [(String,String)]
--this is an inefficient way to check that this alignment will be valid
parseFasta' ((name,""):(name2,seq2):(name3,seq3):xs) bs | (length seq2) /= (length seq3) = error $ "lengths of sequences " ++ name2 ++ " and " ++ name3 ++ " don't match"
parseFasta' old [] =  old
parseFasta' old bs =  case L.unpack (L.take 1 (head bs)) of 
                      ['>'] -> parseFasta' ((trim $ L.unpack (L.drop 1 (head bs)),"") : old) $ tail bs
                      _ -> parseFasta' (appendString old $ trim (L.unpack (head bs))) $ tail bs 

parseFasta :: Monad m => [L.ByteString] -> m [(String,String)]
parseFasta input = return $ parseFasta' [] input 

trim :: String -> String
trim = reverse . dropWhile isSpace . reverse . dropWhile isSpace 

appendString :: [(String,String)] -> String -> [(String, String)]
appendString old add = case  old of 
                (name,x):xs -> (name,x++add):xs

 
