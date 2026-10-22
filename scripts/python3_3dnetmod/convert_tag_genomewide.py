def convert_tag_genomewide(tag, settings):
        pieces = tag.split('_')
        new_tag = ''
        for i in range(len(pieces)-1):
                piece = pieces[i]
                if piece.startswith('chr'):
                        newpiece = 'genomewide'
                else:
                        newpiece = piece
                new_tag= new_tag + newpiece + '_'
        new_tag = new_tag + pieces[-1]
        return new_tag
