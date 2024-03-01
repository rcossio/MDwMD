Next steps:

Improvements:
- Make a table with information from PMID. It is static information and few entries, and making the api request is slow (because of the frequency limit)

Bugfixes:
- If it is loading all the data and they I filter it keep looking for the old entries.. and mixing them with the filtered entries (chatGPT proposed using AbortController)
- results are sorted by accesion number, but the fetch may take longer for some of them and in the end they dont look ordered in the table