import System.Environment (getArgs)
import Phylo.Tree
import System.IO
import Data.Char
import System.Exit

main = do (fn:args) <- getArgs
          file <- openFile fn ReadMode
          contents' <- hGetContents file
          let contents = trim contents'
          case (readLNewickTree contents) of 
                Left err -> do hPutStr stderr err
                               exitWith $ ExitFailure 1
                Right t -> putStrLn $ toNewickL $ unrootL $ prune args t
          return True

trim = f . f where
         f = reverse . dropWhile isSpace
